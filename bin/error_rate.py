#!/usr/bin/env python3

from AmpliconPE import mismatcher
from scipy.stats import trim_mean
import pandas as pd
import numpy as np

DTYPES = dict(
    alignment="O",
    sample="O",
    max_value=np.int64,
    barcode="O",
    reads=np.int64,
    score=np.int64,
)


pileups = (
    pd.read_csv("aligned.csv.gz", dtype=DTYPES, index_col=["sample", "barcode"])[
        "reads"
    ]
    .unstack()
    .T.sort_index()
)

# from AmpliconPE import BarcodeSet
# neighbors = BarcodeSet(pd.Series(np.arange(len(pileups)) ,index=pileups.index, name='iloc'), n_mismatches=1, robust=True, InDels=True)


topN = 100
seed_barcodes = pileups.groupby(level="sample").nlargest(topN)
max_background_barcode_size = 10
null_barcodes = (
    pileups[pileups <= max_background_barcode_size]
    .groupby(level="sample")
    .nlargest(topN)
)


def neighbors(piles):
    return pd.concat(
        {
            barcode: piles.droplevel("sample").reindex(
                set(mismatcher(barcode, 1, InDels=True)) - {barcode}
            )
            for barcode in piles.index.get_level_values("barcode").values
        },
        names=["seed barcode"],
    )


neighboring_barcodes = null_barcodes.groupby(level="sample").apply(neighbors)

background_error_rates = (neighboring_barcodes != 0).groupby(level="sample").mean()

print(f"{len(null_barcodes)} null barcodes have the following background error_rates:")
print(background_error_rates)

assert False

background_rate = background_error_rates.mean()


def trunc_mean(S):
    return trim_mean(S, background_rate)


def trunc_sqrt_mean(S):
    return trim_mean(np.sqrt(S), background_rate) ** 2


def sqrt_mean(S):
    return np.sqrt(S).mean() ** 2


spawns = seed_barcodes.groupby(level="sample").apply(neighbors)

spawn_stats = spawns.groupby(level=["sample", "seed barcode"]).agg(
    [np.median, np.mean, trunc_mean, trunc_sqrt_mean, sqrt_mean],
)

spawn_stats["seed size"] = seed_barcodes

correlations = spawn_stats.groupby(level="sample").corr()["seed size"]
mean_correlations = correlations.unstack().mean()
mean_correlations.pop("seed size")
print(mean_correlations)

nucleotides = list("ACGTN")


def barcode_content(barcodes):
    return pd.concat(
        {nuc: barcodes.str.count(nuc) for nuc in nucleotides}, names=["content"]
    )


spawn_df = spawns.reset_index(["seed barcode", "barcode"])
delta_content = spawn_df["barcode"].apply(barcode_content) - spawn_df[
    "seed barcode"
].apply(barcode_content)

spawns.set_index(delta_content, append=True, inplace=True)
final_ratios = (
    spawns["reads"]
    / seed_barcodes.reindex(spawns["seed barcode"], level="barcode").values
)
final_ratios.name = "Probability"

error_rates = (
    final_ratios.groupby(level=["sample"] + nucleotides)
    .agg(np.mean)
    .reset_index(nucleotides)
)
delta_size = error_rates[nucleotides].sum(axis=1)
deletions = delta_size == -1
insertions = delta_size == +1

assert (error_rates.loc[~(deletions | insertions), nucleotides] == 1).sum(axis=1) == 1
assert (error_rates.loc[~(deletions | insertions), nucleotides] == -1).sum(axis=1) == 1

error_rates["From"] = error_rates[nucleotides].idxmin(axis=1)
error_rates["To"] = error_rates[nucleotides].idxmax(axis=1)

inv_deletion = "Deletion (inverted From/To)"  # Invert Deletions for graphing
error_rates.loc[deletions, "To"] = error_rates.loc[deletions, "From"]
error_rates.loc[deletions, "From"] = inv_deletion
error_rates.loc[insertions, "From"] = "Insertion"

error_rates.to_csv("substitution_error_rates.csv")

import seaborn as sns
import matplotlib.pyplot as plt

ax = plt.gca()
sns.barplot(
    data=error_rates.reset_index(),
    x="from",
    y="Probability",
    hue="to",
    order=nucleotides + ["Insertion" + inv_deletion],
    ax=ax,
)
plt.savefig("Substitution_Error_Rates.pdf")
