#!/usr/bin/env python3
from AmpliconPE import BarcodeSet
from AmpliconPE.shared import DTYPES
from pathlib import Path
import pandas as pd

from jsonargparse import CLI


def error_rate(
    FASTQ_directory: Path,
    align_file: Path = "aligned.csv",
    barcodes_queried: int = 1000,
    max_background_barcode_size: int = 10,
):

    pileups = pd.read_csv(
        FASTQ_directory / align_file, dtype=DTYPES, index_col="barcode"
    )["reads"]
    background_barcodes = pileups[pileups <= max_background_barcode_size].nlargest(
        barcodes_queried
    )

    background_set = frozenset(BarcodeSet(background_barcodes, robust=True).keys())
    all_barcodes = frozenset(pileups.index)
    background_neighbors = len(background_set & all_barcodes)

    largest_barcodes = pileups.nlargest(barcodes_queried)
    largest_set = BarcodeSet(
        zip(largest_barcodes.index, largest_barcodes.index), robust=True
    )

    neighbors = (
        pileups[frozenset(largest_set.keys()) & all_barcodes]
        .reset_index()
        .rename(columns={"reads": "spawn reads", "barcode": "spawn barcode"})
    )
    neighbors["seed barcode"] = neighbors["spawn barcode"].map(largest_set)
    neighbors["spawn reads"] = pileups[neighbors["seed barcode"]]
    neighbors.append(
        {
            "spawn barcode": "Background Tally",
            "spawn reads": len(background_neighbors),
            "seed barcode": "(not read #)",
            "seed reads": len(background_barcodes),
        },
        ignore_index=True,
    ).to_csv("error_rates.csv", index=False)


if __name__ == "__main__":
    CLI()
