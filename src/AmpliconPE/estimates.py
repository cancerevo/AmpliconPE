"""
Functions to estimate technical errors in Amplicon Libraries. 

iHop_rate ~ index-hopping rate between samples
spawn_rate ~ post-transduction SNV rate (this is lower than the SNV rate of barcode flanking sequences 
                                        b/c that error rate includes oligosynthesis+cloning SNVs).
"""

import numpy as np
import pandas as pd


def corr_mean(M):
    """Mean of correlation matrix (off-diagonal)."""
    return np.nanmean(M.values[np.triu_indices(len(M), 1)])


def iHop_rate(barcode_matrix, topN=10):
    """Use the largest barcode pileups and their recurrence in other samples to estimate i-Hopping rate.

    Calculates ratio of largest barcode pileups in each sample to their size in other samples. Computes
    mean of means of these ratios (median of ratios in other samples

    Parameters:
    -----------

    topN (int, default: 10) : Top N barcode pileups from each sample to use."""
    clean = barcode_matrix.loc[
        barcode_matrix.index.get_level_values("target") != "Spike"
    ].T
    largest = clean.max()
    ratios = clean.div(clean).fillna(0)
    S = ratios.T.set_index([clean.idxmax(), largest], append=True).stack()
    S.index.names = ["target", "barcode", "largest_ix", "largest", "Sample"]

    def summary_stats(df, sigmas=5):
        data = df.nlargest(topN, "largest")
        var_stabilized = np.sqrt(data[0] + 0.5)
        no_outliers = data[0].loc[
            var_stabilized < var_stabilized.mean() + sigmas * var_stabilized.std()
        ]
        return (
            no_outliers.apply(["mean", "std"])
            if len(no_outliers) > 1
            else pd.Series(dict(mean=np.nan, std=np.nan))
        )

    summary = (
        S.reset_index("largest")
        .groupby(level=["largest_ix", "Sample"])
        .apply(summary_stats)
    )
    means = summary["mean"]
    M = means.unstack().T
    failure_rate = (means == 0).mean()
    if failure_rate > 0.1:
        print(
            "iHop rate is very low -- cannot estimate for {:.0%} of sample pairs.".format(
                failure_rate
            )
        )
    print(
        """Pearson Correlations:
  iHopping rate between samples: {:.0%}
  Mean Fano Factor: {:.0%}""".format(
            corr_mean(M.corr()), summary.eval("std*std/mean").dropna().mean()
        )
    )
    return means.median()


def spawn_rate(S, frac=0.25, maxN=1000, trim=0.025):
    """Tallies number of read errors from largest tumors to estimate post-transduction error rate.

    Parameters:
    -----------
    frac (float, default: 0.25) : fraction of possible Single Nucleotide Variants (SNVs) that must be observed
                                  (tallies > 0) to warrant use in spawn rate.
    maxN (int, default: 1000)   : Maximum # of tumors to use in estimate.
    trim (float, default: 0.05) : Use trimmed-mean to estimate spawn rate."""
    S = S.loc[S.index.values[0][0]]
    seeds = S.drop("Spike", level="target").nlargest(maxN).sort_values(ascending=False)
    observed = []
    ratios = []
    max_spawn = 3 * len(S.index.values[0][1])
    for i, ((target, barcode), abundance) in enumerate(seeds.iteritems()):
        spawn = S.loc[
            target,
            {
                barcode[:i] + nuc + barcode[i + 1 :]
                for nuc in "ACGT"
                for i, bc_nuc in enumerate(barcode)
                if nuc != bc_nuc
            },
        ]
        if len(spawn) <= frac * max_spawn:
            i -= 1
            break
        observed += spawn.tolist()
        ratios += (spawn / abundance).tolist()

    df = pd.DataFrame(dict(observed=observed, ratios=ratios))

    reduced = df.nlargest(int(i * max_spawn * (1 - trim)), "ratios").nsmallest(
        int(i * max_spawn * (1 - 2 * trim)), "ratios"
    )
    return df["observed"].sum(), seeds.iloc[0:i].sum()
