#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from AmpliconPE.shared import nucleotides, mutation_categories
from pathlib import Path

directory = Path.cwd()
imput_filename = "error_rate.csv.gz"
info_filename = "error_info.csv.gz"
model_filename = "error_model.csv.gz"
truncate = True
transformations = {
    "linear": lambda S: S.mean(),
    "poisson": lambda S: np.mean(np.sqrt(S)) ** 2,
    "log": lambda S: np.exp(np.mean(np.log(S))),
}

transformation = "poisson"
estimator = transformations[transformation]


figure_name = "Substitution_Error_Rates.pdf"
background_info = pd.read_csv(
    directory / info_filename, index_col=["sample", "ensemble"]
)["background"].unstack("ensemble")
background_rates = background_info.eval("`spawn found` / `spawn queried`")


def truncate_sample(df):
    # Add zeros
    sample = df.pop("sample").values[0]
    return df.nsmallest(
        int(round(len(df) * (1 - background_rates[sample]))), "To reads"
    )


error_data = pd.read_csv(directory / imput_filename)

if truncate:
    untruncated = dict(reads=error_data["To reads"].sum(), barcodes=len(error_data))
    error_data = error_data.groupby("sample").apply(truncate_sample)
    truncated = dict(reads=error_data["To reads"].sum(), barcodes=len(error_data))
    ratios = {
        name: 1 - truncated[name] / untruncated[name] for name in untruncated.keys()
    }
    print(
        f"Truncation removed {ratios['reads']:.0%} of reads and {ratios['barcodes']:.0%} of barcodes."
    )


error_rates = (
    error_data.groupby(["sample", "From", "To"])[["From reads", "To reads"]]
    .agg(estimator)
    .eval("Probability = `To reads` / `From reads`")
    .reset_index("sample")
)
error_rates["Category"] = mutation_categories.loc[error_rates.index]

error_rates.reset_index().groupby(["From", "To", "Category"])[
    "Probability"
].mean().to_csv(directory / model_filename)

print("Mean Rates:")
categorical_rates = error_rates.groupby("Category").mean()
categorical_rates["typical"] = pd.Series(
    {
        "read mismatch": 1e-4,
        "transition": 1e-5,
        "transversion": 1e-6,
        "Deletion": 1e-7,
        "Insertion": 1e-8,
    }
)
print(
    categorical_rates.sort_values("typical", ascending=False).to_string(
        float_format="{:.1e}".format
    )
)

sns.set_style("ticks")
rows = 3
fig, axs = plt.subplots(rows, figsize=(9, 6 * rows))

ACGT = list(nucleotides.replace("N", ""))
substitutions = (
    error_rates.groupby(["From", "To"]).mean()["Probability"].unstack().loc[ACGT, ACGT]
)
sns.heatmap(substitutions, square=True, ax=axs[0], cmap="YlOrBr", vmin=0)

sns.barplot(
    data=error_rates.query('Category == "read mismatch"').reset_index(),
    x="From",
    y="Probability",
    order=ACGT,
    ax=axs[1],
)

axs[0].set(title="Cloning+Experiment+PCR Error Rate\n(i.e. Fwd-Rev matching errors)")
axs[1].set(title="Sequencing Error Rate\n(i.e. $to$ `N` base)")
axs[2].set(title="InDels")


error_rates.reset_index(inplace=True)
error_rates["InDel Change"] = error_rates.apply(
    lambda row: row.From
    if row.Category == "Deletion"
    else (row.To if row.Category == "Insertion" else np.nan),
    axis=1,
)

sns.barplot(
    data=error_rates.dropna(),
    x="Category",
    y="Probability",
    hue="InDel Change",
    hue_order=ACGT,
    ax=axs[2],
)
plt.savefig(figure_name)
