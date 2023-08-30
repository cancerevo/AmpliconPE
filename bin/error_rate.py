#!/usr/bin/env python3

from AmpliconPE import mismatcher, barcode_content
from AmpliconPE.shared import DTYPES
import pandas as pd
import numpy as np
from pmap import pmap
from scipy.stats import trim_mean

pileups = pd.read_csv("aligned.csv.gz", dtype=DTYPES, index_col=["sample", "barcode"])[
    "reads"
]


topN = 100
seed_barcodes = pileups.groupby(level="sample").nlargest(topN).droplevel(0)
max_background_barcode_size = 10
null_barcodes = (
    pileups[pileups <= max_background_barcode_size]
    .groupby(level="sample")
    .nlargest(topN)
).droplevel(0)
all_barcodes = pileups.index.get_level_values("barcode").unique()


def get_intersecting_barcodes(barcode):
    return all_barcodes.intersection(
        list(mismatcher(barcode, mismatches=1, InDels=True))
    )


neighboring_barcodes = dict(
    zip(all_barcodes, pmap(get_intersecting_barcodes, all_barcodes))
)


def find_neighbors(df):
    sample = df.index.get_level_values("sample").values[0]
    possible_spawns = pileups.loc[sample]
    all_neighbors = pd.concat(
        {
            barcode: possible_spawns.reindex(neighboring_barcodes[barcode])
            for barcode in df.index.get_level_values("barcode")
        }
    )
    all_neighbors.index.names = ["seed barcode", "spawn barcode"]
    return all_neighbors


neighbors = (
    null_barcodes.groupby(level="sample")
    .apply(find_neighbors)
    .groupby(level=["sample", "seed barcode"])
    .size()
)

background_error_rates = neighbors.groupby(level="sample")["neighbors"].mean()

print(f"{len(null_barcodes)} null barcodes have the following background error_rates:")
print(background_error_rates)
example_bc = pileups.idxmax()[1]
mismatches = list(mismatcher(example_bc, mismatches=1, InDels=True))
print(f"Neighbors per Barcode: {len(mismatches)}")


background_rate = background_error_rates.mean() / len(mismatches)


def trunc_mean(S):
    return trim_mean(S, background_rate)


def trunc_sqrt_mean(S):
    return trim_mean(np.sqrt(S), background_rate) ** 2


def sqrt_mean(S):
    return np.sqrt(S).mean() ** 2


spawns = seed_barcodes.groupby(level="sample").apply(find_neighbors)

spawn_stats = spawns.groupby(level=["sample", "seed barcode"]).agg(
    [np.median, np.mean, trunc_mean, trunc_sqrt_mean, sqrt_mean]
)

spawn_stats["seed reads"] = seed_barcodes

correlations = spawn_stats.groupby(level="sample").corr()["seed reads"]
mean_correlations = correlations.unstack().mean()
mean_correlations.pop("seed reads")
print(mean_correlations)
#
# median             0.467
# mean              -0.229
# trunc_mean         0.775
# trunc_sqrt_mean    0.779
# sqrt_mean         -0.024

nucleotides = list("ACGTN")


seed_barcodes.index.names = ["sample", "seed barcode"]
ratios = (
    spawns.reset_index("spawn barcode")
    .join(
        seed_barcodes,
        lsuffix="_spawn",
        rsuffix="_seed",
    )
    .eval("Probability = reads_spawn / reads_seed")
    .reset_index()
)


content = ratios[["seed barcode", "spawn barcode"]].apply(barcode_content)
delta_content = (content["spawn barcode"] - content["seed barcode"]).unstack("content")


error_rates = (
    ratios.join(delta_content)
    .set_index(nucleotides + ["sample"], append=True)
    .groupby(level=["sample"] + nucleotides)["Probability"]
    .agg(trunc_sqrt_mean)
).reset_index()

error_rates["Type"] = (
    error_rates.reset_index(nucleotides)[nucleotides]
    .sum(axis=1)
    .map({-1: "Deletion", +1: "Insertion", 0: "Substitution"})
)
assert not error_rates["Type"].isnull().any()


error_rates["From"] = error_rates[nucleotides].idxmin(axis=1)
error_rates["To"] = error_rates[nucleotides].idxmax(axis=1)

error_rates.to_csv("substitution_error_rates.csv", index=False)
