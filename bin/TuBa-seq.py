#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, get_PE_FASTQs
import numpy as np
from collections import Counter


FLANK_LENGTH = 8
full_amplicon = "".join(
    (
        "GCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA",  # forward flank
        "........",  # sgID
        "AA.....TT.....AA.....",  # random barcode
        "ATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCT",
    )  # aft flank
)
default_master_read = full_amplicon[
    full_amplicon.find(".")
    - FLANK_LENGTH : full_amplicon.rindex(".")
    + 1
    + FLANK_LENGTH
].replace(".", "N")


def derep(
    FASTQ_directory: Path,
    sgRNA_file: Path = Path("sgRNA_info.csv"),
    master_read: str = default_master_read,
    min_align_score: float = 0.75,
    mismatches_tolerated: int = 1,
):
    """Extracts TuBa-seq double barcodes from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        FASTQ_Directory (or file) containing the forward FASTQ run(s)
        sgRNA_file: All sgRNAs used and their corresponding genes
        master_read: Flanking amplicon sequence to align to each read
        min_align_score: Combined PE alignment score needed to use read, Range [0, 1)
        mismatches_tolerated: # of mismatches tolerated in sgID"""

    sg_info = pd.read_csv(sgRNA_file, converters={"ID": str.upper})
    sgID_length = int(sg_info["ID"].str.len().median())

    sgRNA_map = BarcodeSet(
        sg_info.set_index("ID")["target"],
        n_mismatches=mismatches_tolerated,
        robust=True,
    )

    sgRNA_map.pop_nonunique()

    master_read = MasterRead(master_read)

    pileups = Counter()
    scores = np.zeros(master_read.max_score + 1, dtype=np.int64)
    min_int_score = int(min_align_score * master_read.max_score)

    file_pair = get_PE_FASTQs(FASTQ_directory)
    FASTQ_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in FASTQ_iter:

        double_alignment = master_read.align(fwd_dna, rev_dna)
        scores[double_alignment.score] += 1
        if double_alignment.score <= min_int_score:
            continue

        barcode = double_alignment.extract_barcode()
        if barcode == "Length Mismatch":
            continue

        known_barcode = barcode[:sgID_length]
        sgRNA_target = sgRNA_map.get(known_barcode, "Unknown Target")
        random_barcode = barcode[sgID_length:]

        pileups[(sgRNA_target, random_barcode)] += 1

    poor_alignment = scores[:min_int_score].sum()

    pileups = pd.Series(pileups, name="reads")
    pileups.index.names = "target", "barcode"
    pileups.to_csv(FASTQ_directory / "pileups.csv")

    stats = pd.concat(
        {
            "alignment scores": pd.Series(scores),
            "lost reads": pd.Series(
                {
                    "Poor Alignment": poor_alignment,
                    "Length Mismatch": pileups.sum() - poor_alignment,
                    "Index Mismatches": FASTQ_iter.index_mismatches,
                }
            ),
        }
    )
    stats.index.names = "type", "outcome"
    stats.name = "reads"
    stats.to_csv(FASTQ_directory / "stats.csv")


if __name__ == "__main__":
    CLI()
