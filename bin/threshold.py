#!/usr/bin/env python3
"""
"""
import pandas as pd
from pathlib import Path
from jsonargparse import CLI
from AmpliconPE.shared import DTYPES

threshold = 0.5


def threshold(
    directory: Path = "consolidated",
    threshold: float = 0.5,
):

    info = (
        pd.read_csv(
            directory / "info.csv.gz",
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

    pileups = pd.read_csv(directory / "pileups.csv.gz", dtype=DTYPES)
    gb = pileups.set_index("score").groupby(
        lambda ix: "aligned" if ix > scaled_threshold else "unaligned"
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
                "aligned": total_reads["aligned"],
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
{fractions['aligned']:.1%} of reads align to the amplicon when a {threshold:.1%} normalized score threshold is applied.
"""
    )

    outcomes = ["Adapter Index Mismatches", "unaligned", "aligned"]

    if index_mismatches.sum() == 0:
        print(f"No {outcomes[0]}.")
        outcomes.pop(0)
    else:
        df = sample_reads.unstack(level="outcome")
        df[outcomes[0]] = index_mismatches

    print(sample_reads[outcomes].to_string())


if __name__ == "__main__":
    CLI()
