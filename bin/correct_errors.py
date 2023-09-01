import pandas as pd
from AmpliconPE.shared import barcode_content, mismatcher

nucleotides = list("ACGTN")
error_rate_file = "error_rates.csv"

error_rates = pd.read_csv(error_rate_file)
simplified = error_rates.groupby("Type")["Probability"].mean()

gb = error_rates.groupby(nucleotides)[["From", "To"]]
content_mapper = gb.value_counts()
assert len(gb) == len(content_mapper)

aligned_file = "aligned.csv"
aligned = pd.read_csv(aligned_file, index_col="barcode")
barcodes = aligned.index

content = barcodes.apply(barcode_content)
nucletoide_change = content.apply(content_mapper)


def simplified_rate(df):
    return simplified[df["Type"].iloc[0]]


error_model = error_rates.groupby(nucleotides).agg(simplified_rate)


def spawn_reads(barcode):
    neighbors = barcodes.intersection(
        list(mismatcher(barcode, mismatches=1, InDels=True))
    )
    barcode_content = content[barcode]
    delta_content = content[neighbors] - barcode_content
    return (error_model[delta_content] * aligned.loc[neighbors, "reads"]).sum()


aligned["Spawn Reads"] = barcodes.apply(spawn_reads)
aligned["Seed Reads"] = (error_model * content).sum(axis=1) * aligned["reads"]

aligned.to_csv("error_corrected.csv")
