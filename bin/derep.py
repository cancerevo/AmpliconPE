#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import SimplexMasterRead as MasterRead, pairedFASTQiter, get_PE_FASTQs
from collections import Counter


def derep(FASTQ_directory: Path, master_read: str, trim_fraction: float = 0.75):
    """Lossless barcode extraction from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        FASTQ_directory (or file) containing the forward/reverse-read FASTQs
        master_read: Flanking amplicon sequence to align to each read
        trim_fraction: Fraction of barcode flanks to trim before core alignment."""

    # in case running interactively, convert strings to paths
    FASTQ_directory = Path(FASTQ_directory)

    master_read = MasterRead(master_read, trim_fraction=trim_fraction)

    pileups = Counter()
    alignment_scores = Counter()

    file_pair = get_PE_FASTQs(FASTQ_directory)
    FASTQ_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in FASTQ_iter:

        alignment = master_read.align(fwd_dna, rev_dna)
        alignment_scores[alignment.get_scores()] += 1
        pileups[(alignment.final_score, alignment.extract_barcode())] += 1

    ## Output

    self_alignment = master_read.self_alignment
    align_pairs = self_alignment.align_pairs

    info = dict(
        zip(
            map("{:} max score".format, align_pairs),
            master_read.self_alignment.get_scores(),
        )
    )
    info["index mismatches"] = FASTQ_iter.index_mismatches

    index_labels = dict(
        pileups=["score", "barcode"], alignment_scores=align_pairs, info=["stat"]
    )

    for name, index_label in index_labels.items():
        pd.Series(eval(name), name="reads").to_csv(
            FASTQ_directory / f"{name}.csv", index_label=index_labels[name]
        )


if __name__ == "__main__":
    CLI()
