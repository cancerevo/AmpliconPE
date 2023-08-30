#!/usr/bin/env python3
"""
"""
import pandas as pd
from pathlib import Path
from jsonargparse import CLI
from AmpliconPE.shared import DTYPES


def threshold(
    directory: Path,
    pileup_file: Path = "pileups.csv",
    info_file: Path = "info.csv",
    threshold: float = 0.5,
    output_names: list = ["unaligned", "aligned"],
):

    info = pd.read_csv(directory / info_file, dtype=DTYPES).set_index("alignment")
    max_score = info.loc["core-ref", "max_value"]

    scaled_threshold = int(max_score * threshold)

    pileups = pd.read_csv(directory / pileup_file, dtype=DTYPES).eval(
        f"""
        passed_threshold = score > {scaled_threshold}
        total_score = score*reads
        perfect = (score == max_score)*reads"""
    )

    totals = pileups.groupby(["passed_threshold", "barcode"])[
        ["total_score", "reads", "perfect"]
    ].sum()
    totals["mean alignment"] = totals.pop("total_score") / (totals.reads * max_score)
    info["perfect reads"] = totals.pop("perfect").sum()

    for threshold_status, df in pileups.groupby(level="aligned"):
        name = output_names[int(threshold_status)]
        df.to_csv(f"{name}.csv")
        info[f"{name} reads"] = df.reads.sum()

    info.to_csv("info.csv")


if __name__ == "__main__":
    CLI()
