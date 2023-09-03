#!/usr/bin/env python3
from AmpliconPE import BarcodeSet
from AmpliconPE.shared import DTYPES, annotate_change, barcode_content
import pandas as pd
import numpy as np
from scipy.stats import trim_mean

Directory = ""
aligned_file = "aligned.csv.gz"
pileups = pd.read_csv(Directory, dtype=DTYPES, index_col="barcode")["reads"]

topN = 100
max_background_barcode_size = 10

background_barcodes = pileups[pileups <= max_background_barcode_size].nlargest(topN)

background_set = frozenset(BarcodeSet(background_barcodes, robust=True).keys())
all_barcodes = frozenset(pileups.index)
neighbors = len(background_set & all_barcodes)

background_rate = len(neighbors) / len(background_barcodes)


frac_background_rate = len(neighbors) / len(background_set)


def trunc_sqrt_mean(S):
    return trim_mean(np.sqrt(S), frac_background_rate) ** 2


largest_barcodes = pileups.nlargest(topN)
largest_set = BarcodeSet(
    zip(largest_barcodes.index, largest_barcodes.index), robust=True
)

largest = (
    pileups[frozenset(largest_set.keys()) & all_barcodes]
    .reset_index()
    .rename(columns={"reads": "spawn reads", "barcode": "spawn barcode"})
)
largest["seed barcode"] = largest["spawn barcode"].map(largest_set)
largest["Probability"] = largest["spawn reads"] / pileups[largest["seed barcode"]]

delta_content = (
    barcode_content(largest["spawn barcode"]) - barcode_content(largest["seed barcode"])
).unstack("content")


# Should AGG by K80 model here

#    largest.join(delta_content)
#    .groupby(nucleotides)["Probability"]
error_rates = (
    largest.groupby(delta_content)["Probability"].agg(trunc_sqrt_mean)
).reset_index()

nucleotides = delta_content.columns
annotated_rates = error_rates.join(
    annotate_change(error_rates.reset_index()[nucleotides])
)
annotated_rates.to_csv("error_rates.csv", index=False)

summary = annotated_rates.groupby(level="Type").mean()
print("Overall Error Rates:")
print(summary.to_string())
