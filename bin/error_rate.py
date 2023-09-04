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
    output_file: Path = "error_rate.csv",
    info_file: Path = "error_info.csv",
):

    pileups = pd.read_csv(
        FASTQ_directory / align_file, dtype=DTYPES, index_col="barcode"
    )["reads"]
    all_barcodes = frozenset(pileups.index)

    barcode_sets = dict(
        background=pileups[pileups <= max_background_barcode_size], largest=pileups
    )

    barcode_groups = {}
    for group, barcode_set in barcode_sets.items():
        queried = barcode_set.nlargest(barcodes_queried)
        possible_neighbors = BarcodeSet(zip(queried.index, queried.index), robust=True)
        possible_neighbors.pop_nonunique()
        neighbors = frozenset(possible_neighbors.keys()) & all_barcodes
        barcode_groups[group] = {
            "seeds queried": frozenset(possible_neighbors.values()),
            "seed->spawn": possible_neighbors,
            "spawn found": neighbors,
        }

    largest = barcode_groups["largest"]
    error_df = (
        pileups[largest["spawn found"]]
        .reset_index()
        .rename(columns={"reads": "spawn reads", "barcode": "spawn barcode"})
    )
    error_df["seed barcode"] = error_df["spawn barcode"].map(largest["seed->spawn"])
    error_df["spawn reads"] = pileups[error_df["seed barcode"]]
    error_df.to_csv(FASTQ_directory / output_file, index=False)

    info = pd.DataFrame(
        {name: len(group) for name, group in barcode_groups.items()}
    ).rename({"seed mapping": "spawn queried"}, axis=0)
    info.index.names = ["# of barcodes"]
    info.to_csv(FASTQ_directory / info_file)


if __name__ == "__main__":
    CLI()
