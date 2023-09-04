#!/usr/bin/env python3

from AmpliconPE import BarcodeSet, barcode_content
from AmpliconPE.shared import DTYPES, nucleotides
import pandas as pd
import numpy as np
from pmap import pmap
from scipy.stats import trim_mean

directory = "consolidated"
aligned_file = "aligned.csv.gz"

error_info_file = "error_info.csv"
error_rate_file = "error_rate.csv"
truncate = True
info = pd.read_csv(directory / error_info_file, index_col=["sample", "# of barcodes"])
background_rates = (
    info["background"].groupby(level="sample")["spawn found"].sum()
)  # .agg(lambda S: S['spawn found'] / S['seeds queried'])


def truncate_sample(df):
    # Add zeros
    sample = df["sample"].values[0]
    return df.nsmallest(int(round(len(df) - background_rates[sample])), "spawn reads")


error_df = pd.read_csv(directory / error_rate_file)

if truncate:
    error_df = (
        error_df.groupby(["sample", "seed barcode"])
        .apply(truncate_sample)
        .reset_index()
    )


" BARCODE SET COULD have annotations of every mut built-in. "


def trunc_mean(S):
    return trim_mean(S, frac_background_rate)


def nsmallest(S):
    return S.nsmallest(len(S) - int(round(background_rate))).mean()


gb = spawn.groupby("seed barcode")
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


annotated_rates = error_rates.join(
    annotate_change(error_rates.reset_index(nucleotides)[nucleotides])
).set_index(["Type", "From", "To"] + nucleotides)
annotated_rates.to_csv("error_rates.csv")

summary = annotated_rates.groupby(level="Type").mean()
print("Overall Error Rates:")
print(summary.to_string())
