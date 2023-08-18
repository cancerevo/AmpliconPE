#!/usr/bin/env python3
"""
"""
import pandas as pd
import numpy as np

threshold = 0.5
DTYPES = dict(
    alignment="O",
    sample="O",
    max_value=np.int64,
    barcode="O",
    reads=np.int64,
    score=np.int64,
)


info = (
    pd.read_csv(
        "consolidated_info.csv.gz",
        index_col=["alignment", "sample"],
        usecols=["alignment", "sample", "max_value"],
        dtype=DTYPES,
    )["max_value"]
    .unstack()
    .T
)
index_mismatches = info.pop("Index Mismatches")

max_scores = info["core-ref"].value_counts()
assert len(max_scores) == 1, "Samples have different maximum Core-ref values."
max_score = max_scores.index.values[0]

scaled_threshold = int(max_score * threshold)

pileups = pd.read_csv("consolidated_pileups.csv.gz", dtype=DTYPES)
gb = pileups.set_index("score").groupby(
    lambda ix: "Aligned" if ix > scaled_threshold else "Unaligned"
)

outputs = {}
for quality, pileup_set in gb:
    S = pileup_set.groupby(["sample", "barcode"]).apply(
        lambda df: pd.Series(
            {
                "mean alignment": df.eval("score*reads").sum()
                / (df.reads.sum() * max_score),
                "reads": df.reads.sum(),
            }
        )
    )
    S.to_csv(f"{quality}.csv.gz")
    outputs[quality] = S

sample_reads = (
    pd.concat(outputs, names=["outcome"])
    .groupby(level=["outcome", "sample"])["reads"]
    .sum()
)
total_reads = sample_reads.groupby(level="outcome").sum()

fractions = (
    pd.Series(
        {
            "Aligned": total_reads["Aligned"],
            "perfect": pileups.query("score == @max_score")["reads"].sum(),
        }
    )
    / total_reads.sum()
)

print(
    f"""
Summary Stats:
-------------

{fractions['perfect']:.1%} of reads perfectly match the Master Read.
{fractions['Aligned']:.1%} of reads align to the amplicon when a {threshold:.1%} normalized score threshold is applied.
"""
)

outcomes = ["Adapter Index Mismatches", "Unaligned", "Aligned"]

if index_mismatches.sum() == 0:
    print(f"No {outcomes[0]}.")
    outcomes.pop(0)
else:
    df = sample_reads.unstack(level="outcome")
    df[outcomes[0]] = index_mismatches


print(sample_reads[outcomes].to_string())

# sns.barplot(
#    data=sample_reads.reset_index(),
#    x="sample",
#    y="reads",
#    hue="outcome",
#    hue_order=outcomes,
#    ax=axs[3],
# )
