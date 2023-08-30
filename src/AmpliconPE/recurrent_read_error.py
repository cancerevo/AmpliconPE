seed_error_model = error_model_inverted_InDels.groupby("From").sum()


def spawn_reads(barcode):
    neighbors = barcodes.intersection(list(mismatcher(mismatches=1, InDels=True)))
    barcode_content = content[barcode]
    delta_content = content[neighbors] - barcode_content
    return (error_model[delta_content] * aligned.loc[neighbors, "reads"]).sum()


corrected_barcodes = aligned.assign(
    {
        "Spawn Reads": lambda df: spawn_reads(df["barcode"]),
        "Seed Reads": lambda df: (seed_error_model * content).sum() * df["reads"],
    }
)

corrected_barcodes.to_csv("error_corrected.csv")
