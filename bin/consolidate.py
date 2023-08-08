#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path

table_index_names = dict(pileups=["target", "barcode"], stats=["type", "outcome"])


def consolidate(
    directory: Path,
    glob: str = "*",
    consolidated_prefix: str = "consolidated_",
    destructive: bool = False,
    ignore_combined: bool = True,
):
    """Globs csv output files into a single csv.

    Expects the following Directory Tree:

    <directory>
     |- <Sample directory>
         |- pileups.csv
         |- <any output>.csv

    Args:
        directory: Directory containing sample directories
        glob: csv basename to glob (e.g. `pileups`, `stats`, ...)
        consolidated_prefix: prefix of basename of consolidated csv(s).
        destructive: Remove original files."""

    if not directory.exists():
        raise FileNotFoundError(f"{directory} does not exist.")
    if not directory.is_dir():
        raise ValueError(f"{directory} is not a directory.")

    globs = list(directory.glob(f"**/{glob}.csv"))
    basenames = set(map(lambda p: p.name, globs))
    samples = set(map(lambda p: p.parts[-2], globs))

    print(f"Found {len(globs)} files.")
    print(f"Found the following {len(basenames)} basenames:")
    print(*basenames)
    print(f"Found the following {len(samples)} sample directories:")
    print(*samples)

    for basename in basenames:
        dfs = {}
        for sample in samples:
            try:
                filename = directory / sample / basename
                dfs[sample] = pd.read_csv(filename)
                if destructive:
                    filename.unlink()
            except FileNotFoundError:
                print(f"Could not find {filename}")

        combined = pd.concat(dfs, names=["sample"])
        combined.reset_index(0).to_csv(
            Path(consolidated_prefix + basename), index=False
        )


if __name__ == "__main__":
    CLI()
