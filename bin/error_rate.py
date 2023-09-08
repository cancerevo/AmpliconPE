#!/usr/bin/env python3
from AmpliconPE import BarcodeSet
from AmpliconPE.shared import DTYPES, annotate_mutations
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

    barcode_group = {}
    for group, barcode_set in barcode_sets.items():
        queried = barcode_set.nlargest(barcodes_queried)
        possible_neighbors = BarcodeSet(
            zip(queried.index, queried.index),
            robust=True,
            InDels=True,
            fixed_length=False,
        )
        usable_neighbors = {
            k: v for k, v in possible_neighbors.items() if k != v and type(v) != set
        }
        barcode_group[group] = usable_neighbors

    all_mutations = (
        pd.Series(barcode_group["largest"], name="From")
        .reset_index()
        .rename(columns=dict(index="To"))
    )
    annotations = annotate_mutations(all_mutations)
    for name, barcodes in all_mutations.iteritems():
        annotations[f"{name} reads"] = barcodes.map(pileups).fillna(0).astype(int)

    annotations.sort_values("From reads").to_csv(
        FASTQ_directory / output_file, index=False
    )

    info = pd.concat(
        {
            name: pd.Series(
                {
                    "seeds queried": len(frozenset(mapping.values())),
                    "spawn queried": len(mapping),
                    "spawn found": len(all_barcodes & frozenset(mapping.keys())),
                },
                name="barcodes",
            )
            for name, mapping in barcode_group.items()
        }
    )
    info.index.names = ["group", "ensemble"]
    info.unstack().T.to_csv(FASTQ_directory / info_file)


if __name__ == "__main__":
    CLI()
