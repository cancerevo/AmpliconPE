#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path

table_index_names = dict(pileups=["target", "barcode"], stats=["type", "outcome"])


def consolidate(
    directory: Path,
    pileups: str = "pileups.csv",
    stats: str = "stats.csv",
    output: Path = Path("consolidated_output.h5"),
    destructive: bool = False,
):
    """Extracts TuBa-seq double barcodes from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        directory: Directory (or file) containing directories of pileups/stats to consolidate.
        pileups: Name of pileup files
        stats: Name of stats files
        output: Name of HDF5 output store
        destructive: Remove pileups.csv & stats.csv"""

    tables = dict(
        pileups=pd.concat(
            {
                filename.parent.name: pd.read_csv(
                    filename, index_col=["target", "barcode"]
                )["reads"]
                for filename in directory.glob(f"**/{pileups}")
            },
            names=["sample"],
        ),
        stats=pd.concat(
            {
                filename.parent.name: pd.read_csv(
                    filename, index_col=["type", "outcome"]
                )["reads"]
                for filename in directory.glob(f"**/{stats}")
            },
            names=["sample"],
            axis=1,
        ),
    )

    with pd.HDFStore(output, "w", complevel=9) as store:
        for name, table in tables.items():
            store.put(name, table)


if __name__ == "__main__":
    CLI()
