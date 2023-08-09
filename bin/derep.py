#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import SimplexMasterRead as MasterRead, pairedFASTQiter, get_PE_FASTQs
from collections import Counter


def derep(FASTQ_directory: Path, master_read: str, trim_fraction: float = 0.6):
    """Lossless barcode extraction from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        FASTQ_directory (or file) containing the forward FASTQ run(s)
        master_read: Flanking amplicon sequence to align to each read
        trim_fraction: Fraction of barcode flanks to trim before core alignment."""

    # in case running interactively, convert strings to paths
    FASTQ_directory = Path(FASTQ_directory)

    master_read = MasterRead(master_read, trim_fraction=trim_fraction)

    pileups = Counter()
    scores = Counter()

    file_pair = get_PE_FASTQs(FASTQ_directory)
    FASTQ_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in FASTQ_iter:

        alignment = master_read.align(fwd_dna, rev_dna)
        scores[alignment.get_scores()] += 1
        pileups[(alignment.final_score, alignment.extract_barcode())] += 1

    ## Output
    max_scores = master_read.self_alignment.get_scores()

    score_series = pd.Series(
        scores.values(),
        index=pd.MultiIndex.from_tuples(
            scores.keys(),
            names=[
                f"{name} (max={max_score})"
                for name, max_score in zip(
                    master_read.self_alignment.align_pairs, max_scores
                )
            ],
        ),
        name="reads",
    )

    score_series.to_csv(FASTQ_directory / "scores.csv")

    pileups = pd.Series(
        pileups.values(),
        name="reads",
        index=pd.MultiIndex.from_tuples(
            pileups.keys(), names=[f"score (max={master_read.max_score})", "barcode"]
        ),
    )
    pileups.to_csv(FASTQ_directory / "pileups.csv")


if __name__ == "__main__":
    CLI()
