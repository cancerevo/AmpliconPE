#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI
from pathlib import Path


def consolidate(
    directory: Path,
    glob: str = "*",
    output_dir: Path = "consolidated",
    destructive: bool = False,
    ignore: list = ["Undetermined", "undetermined"],
    compression: str = ".gz",
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
        output_dir: Directory of consolidated csv(s).
        destructive: Remove original files.
        ignore: Sample names to not consolidate.
        compression:"""

    if not directory.exists():
        raise FileNotFoundError(f"{directory} does not exist.")
    if not directory.is_dir():
        raise ValueError(f"{directory} is not a directory.")
    output_dir.mkdir(exist_ok=True)

    globs = list(directory.glob(f"**/{glob}.csv"))
    basenames = set(map(lambda p: p.name, globs))
    samples = set(map(lambda p: p.parts[-2], globs)) - set(ignore)

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
            output_dir / (basename + compression), index=False
        )


if __name__ == "__main__":
    CLI()
