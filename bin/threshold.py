#!/usr/bin/env python3
"""

Questions This Report Answers:

1. Were my Fwd & Rev reads good quality?
2. Did those two reads cover the same sequence?
3. Is this sequence my barcodes?
4. What else?
    
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
info.pop("Index Mismatches")

max_scores = info["core-ref"].value_counts()
assert len(max_scores) == 1, "Samples have different maximum Core-ref values."
max_score = max_scores.index.values[0]

scaled_threshold = int(max_score * threshold)

pileups = pd.read_csv("consolidated_pileups.csv.gz", dtype=DTYPES)
gb = pileups.set_index("score").groupby(
    lambda ix: "Aligned" if ix > scaled_threshold else "Unaligned"
)

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
