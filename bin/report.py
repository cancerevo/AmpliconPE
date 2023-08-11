#!/usr/bin/env python3
"""

Questions This Report Answers:

1. Were my Fwd & Rev reads good quality?
2. Did those two reads cover the same sequence?
3. Is this sequence my barcodes?
4. What else?
    
"""
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import entropy


threshold = 0.5
HISTOGRAM_BINS = 200

rows = 3
fig, axs = plt.subplots(rows, figsize=(8, 6 * rows))

info = pd.read_csv("consolidated_info.csv", index_col=["alignment", "sample"])[
    "max_value"
]
index_mismatches = info.pop("Index Mismatches")

gb = info.groupby(level="alignment")
max_scores = gb.median().astype(int)
max_final_score = max_scores["core-ref"]
if not np.allclose(gb.std(), np.zeros(len(gb))):
    raise ValueError("Samples have different maximum alignment values.")

align_pairs = max_scores.index.values


pileups = pd.read_csv("consolidated_pileups.csv").eval(
    "total_alignment = score*reads/@max_final_score"
)

pileups_marginal = pileups.groupby(["sample", "barcode"])[
    ["total_alignment", "reads"]
].sum()
pileups_marginal["mean alignment"] = (
    pileups_marginal.pop("total_alignment") / pileups_marginal["reads"]
)


scores = pd.read_csv(
    "consolidated_scores.csv", index_col=["sample"] + list(align_pairs), dtype=np.uint64
)["reads"]


base_scores = {name: scores.groupby(level=name).sum() for name in align_pairs}

normed_scores = {
    name: pd.Series(
        base_score.values,
        index=pd.Index(base_score.index / max_scores[name], name="score"),
        name="reads",
    )
    for name, base_score in base_scores.items()
}


pileup_counts, pileup_bins = np.histogram(
    pileups_marginal["mean alignment"],
    weights=pileups_marginal["reads"],
    bins=HISTOGRAM_BINS,
    range=(0, 1),
)

normed_scores["pileups"] = pd.Series(
    pileup_counts,
    index=pd.Index(pileup_bins[:-1] + pileup_bins[1] / 2, name="score"),
    name="reads",
)
normed_scores = pd.concat(normed_scores, names=["alignment"])

sns.ecdfplot(
    data=normed_scores.reset_index(),
    x="score",
    weights=normed_scores.name,
    hue="alignment",
    hue_order=align_pairs,
    ax=axs[0],
)
axs[0].set(xlim=[0, 1])


def length_stdev(bc):
    return pd.Series(bc).str.len().std()


pile_stats = (
    pileups.groupby("score")
    .agg(dict(barcode=length_stdev, reads=entropy))
    .rename(columns=dict(barcode="Barcode Length StDev", reads="Shannon Entropy"))
)
pile_stats["fwd-rev alignment"] = (
    scores.reset_index().groupby("core-ref")["fwd-rev"].mean() / max_scores["fwd-rev"]
)
pile_stats.columns.name = "stat"
pile_stats /= pile_stats.max() / 100

line_data = pile_stats.groupby(lambda ix: 2 * round(ix / 2)).mean().stack()
line_data.name = "% of Max"

sns.lineplot(
    data=line_data.reset_index(), x="score", y=line_data.name, hue="stat", ax=axs[1]
)

axs[1].axvline(threshold * max_final_score, ls=":")

scaled_threshold = max_final_score * threshold
unaligned = (
    pileups.query("score <= @scaled_threshold")
    .groupby(["sample", "barcode"])[["total_alignment", "reads"]]
    .sum()
)
unaligned["mean alignment"] = unaligned.pop("total_alignment") / unaligned["reads"]
unaligned.sort_values(by="reads", ascending=False).to_csv("unaligned.csv")

aligned = (
    pileups.query("score > @scaled_threshold")
    .groupby(["sample", "barcode"])[["total_alignment", "reads"]]
    .sum()
)
aligned["mean alignment"] = aligned.pop("total_alignment") / aligned["reads"]
aligned.to_csv("aligned_pileups.csv")

aligned_reads = aligned["reads"].sum()
total_reads = aligned_reads + unaligned["reads"].sum()
aligned_frac = aligned_reads / total_reads
perfect_frac = pileups.query("score == @max_final_score")["reads"].sum() / total_reads

print(
    f"""
Summary Stats:
-------------

{perfect_frac:.1%} of reads perfectly match the Master Read.
{aligned_frac:.1%} of reads align to the amplicon when a {threshold:.1%} normalized score threshold is applied.
"""
)


outcomes = ["Adapter Index Mismatches", "Unaligned", "Aligned"]
sample_stats = pd.concat(
    [
        index_mismatches,
        unaligned.groupby(level="sample")["reads"].sum(),
        aligned.groupby(level="sample")["reads"].sum(),
    ],
    keys=outcomes,
    names=["outcome"],
)
sample_stats.name = "reads"

sns.barplot(
    data=sample_stats.reset_index(),
    x="sample",
    y=sample_stats.name,
    hue="outcome",
    hue_order=outcomes,
    ax=axs[2],
)

plt.savefig("report.pdf")
