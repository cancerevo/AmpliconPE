#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

inv_deletion = "Deletion (inverted From/To)"  # Invert Deletions for graphing
error_rates = pd.read_csv(
    "substitution_error_rates.csv", usecols=["sample", "From", "To", "Probability"]
)
sns.set_style("ticks")
rows = 3
fig, axs = plt.subplots(rows, figsize=(9, 6 * rows))

ACGT = list("ACGT")
substitutions = (
    error_rates.groupby(["From", "To"]).mean()["Probability"].unstack().loc[ACGT, ACGT]
)
sns.heatmap(substitutions, square=True, ax=axs[0], cmap="YlOrBr", vmin=0)

sns.barplot(
    data=error_rates.query('To == "N"').reset_index(),
    x="From",
    y="Probability",
    order=ACGT,
    ax=axs[1],
)

axs[0].set(title="Cloning+Experiment+PCR Error Rate\n(i.e. Fwd-Rev matching errors)")
axs[1].set(title="Sequencing Error Rate\n(i.e. $to$ `N` base)")

InDels = ["Insertion", inv_deletion]
sns.barplot(
    data=error_rates.query("From in @InDels").reset_index(),
    x="From",
    y="Probability",
    hue="To",
    hue_order=ACGT,
    ax=axs[2],
)
plt.savefig("Substitution_Error_Rates.pdf")
