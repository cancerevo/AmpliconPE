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

background_error_rates = neighbors.groupby(level="sample").mean()

print(f"{len(null_barcodes)} null barcodes have the following background error_rates:")
print(background_error_rates)
example_bc = pileups.idxmax()[1]
mismatches = list(mismatcher(example_bc, mismatches=1, InDels=True))
print(f"Neighbors per Barcode: {len(mismatches)}")


background_rate = background_error_rates.mean()
frac_background_rate = background_rate / len(mismatches)


def trunc_mean(S):
    return trim_mean(S, frac_background_rate)


def nsmallest(S):
    return S.nsmallest(len(S) - int(round(background_rate))).mean()


spawns = seed_barcodes.groupby(level="sample").apply(find_neighbors)

gb = spawns.groupby(level=["sample", "seed barcode"])
transformations = {
    "none": [lambda S: S, lambda S: S],
    "poisson": [lambda S: np.sqrt(S), lambda S: S**2],
    "log": [lambda S: np.log(S), lambda S: np.exp(S)],
}
spawn_stats = pd.concat(
    {
        (name, stat.__name__): inverse(transformation(gb.agg(stat)))
        for name, (transformation, inverse) in transformations.items()
        for stat in [np.mean, trunc_mean, nsmallest]
    },
    names=["transformation", "estimator"],
).unstack(["transformation", "estimator"])


def correlation(col):
    return pd.concat([col, seed_barcodes]).corr().iloc[0, 1]


correlations = spawn_stats.apply(correlation)
print(correlations)
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


def trunc_sqrt_mean(S):
    return trunc_mean(np.sqrt(S)) ** 2


error_rates = (
    ratios.join(delta_content)
    .set_index(nucleotides + ["sample"], append=True)
    .groupby(level=["sample"] + nucleotides)["Probability"]
    .agg(trunc_sqrt_mean)
).reset_index()


def annotate_change(nucleotide_change):
    K80 = """
* T V V N
T * V V N
V V * T N
V V T * N
N N N N *"""
    K_categories = {
        "T": "transition",
        "V": "transversion",
        "N": "sequencing error",
        "*": "*",
    }
    nucleotides = list("AGCTN")
    mut_categories = (
        pd.DataFrame(
            [list(map(str.split, line)) for line in K80.splitlines()[1:]],
            columns=pd.Index(nucleotides, name="From"),
            index=pd.Index(nucleotides, name="To"),
        )
        .stack()
        .map(K_categories)
    )

    From = nucleotide_change.idxmin(axis=1)
    To = nucleotide_change.idxmax(axis=1)

    Type = nucleotide_change.sum(axis=1).map(
        {-1: "Deletion", +1: "Insertion", 0: "Substitution"}
    )
    out = pd.DataFrame({"From": From, "To": To, "Type": Type})
    out.loc[:, "Substitution"] = out.loc[:, "Substitution"][["From", "To"]].map(
        mut_categories
    )
    assert not out["Type"].isnull().any()
    return out


annotated_rates = error_rates.join(
    annotate_change(error_rates.reset_index(nucleotides)[nucleotides])
).set_index(["Type", "From", "To"] + nucleotides)
annotated_rates.to_csv("error_rates.csv")

summary = annotated_rates.groupby(level="Type").mean()
print("Overall Error Rates:")
print(summary.to_string())
