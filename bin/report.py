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
from AmpliconPE.shared import DTYPES

alignments = ["fwd-ref", "rev-ref", "fwd-rev", "core-ref"]

sns.set_style("ticks")
rows = 3
fig, axs = plt.subplots(rows, figsize=(12, 9 * rows))

info = (
    pd.read_csv(
        "consolidated_info.csv.gz",
        index_col=["stat", "sample"],
        usecols=["stat", "sample", "max_value"],
        dtype=DTYPES,
    )["value"]
    .unstack()
    .T
)

index_mismatches = info.pop("Index Mismatches")

max_scores = info.median().astype(int)

max_final_score = max_scores["core-ref"]
if not np.allclose(info.std(), np.zeros(len(alignments))):
    raise ValueError("Samples have different maximum alignment values.")


pileups = pd.read_csv("consolidated_pileups.csv.gz", dtype=DTYPES).eval(
    "total_alignment = score*reads/@max_final_score"
)

pileups_marginal = pileups.groupby(["sample", "barcode"])[
    ["total_alignment", "reads"]
].sum()
pileups_marginal["mean alignment"] = (
    pileups_marginal.pop("total_alignment") / pileups_marginal["reads"]
)

score_ix = ["sample"] + alignments
scores = pd.read_csv(
    "consolidated_scores.csv.gz",
    index_col=score_ix,
    dtype=DTYPES,
    usecols=score_ix + ["reads"] + alignments,
)["reads"]


base_scores = {name: scores.groupby(level=name).sum() for name in alignments}

normed_scores = {
    name: pd.Series(
        base_score.values,
        index=pd.Index(base_score.index / max_scores[name], name="score"),
        name="reads",
    )
    for name, base_score in base_scores.items()
}

normed_scores = pd.concat(normed_scores, names=["alignment"])

weighted_means = (
    normed_scores.reset_index("score")
    .eval("score_reads = score * reads")[["score_reads", "reads"]]
    .groupby(level="alignment")
    .sum()
    .eval("score_reads / reads")
)

print("Mean Scores:")
print("------------")
print(weighted_means[alignments].to_string(float_format="{:.1%}".format))


sns.histplot(
    data=normed_scores.loc[alignments[:2]].reset_index(),
    x="score",
    weights=normed_scores.name,
    hue="alignment",
    hue_order=alignments[:2],
    ax=axs[0],
    bins=max_scores["fwd-ref"] + 1,
    element="step",
)
axs[0].set(xlim=[0, 1])

##################################
# Self-Ref Corr
##################################

final_scores = scores.groupby(level=alignments[2:]).sum().reset_index()
ax = axs[1]
hb = ax.hexbin(
    final_scores["fwd-rev"].values,
    final_scores["core-ref"].values,
    final_scores["reads"],
    gridsize=final_scores.agg(lambda s: len(set(s)))[alignments[2:]].values,
    cmap="YlOrBr",
    bins="log",
)
ax.set(xlabel="Fwd-Rev Alignment Score", ylabel="Core-Ref Alignment Score")


def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)


def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)


def corr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))


final_corr = corr(
    final_scores["fwd-rev"], final_scores["core-ref"], final_scores["reads"]
)
ax.text(
    0.02,
    0.95,
    f"{final_corr:.1%}",
    horizontalalignment="left",
    verticalalignment="center",
    transform=ax.transAxes,
)

cb = fig.colorbar(hb, ax=ax, label="reads")


def barcode_length_var(df):
    barcode_lengths = df["barcode"].str.len()
    return cov(barcode_lengths, barcode_lengths, df.reads)


bc_length_std = np.sqrt(pileups.groupby("score").apply(barcode_length_var))
bc_length_std.name = "Length St.Dev."

sns.pointplot(
    data=bc_length_std.reset_index(),
    x="score",
    y=bc_length_std.name,
    join=False,
    ax=axs[2],
)


plt.savefig("report.pdf")
