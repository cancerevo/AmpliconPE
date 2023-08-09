#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, get_PE_FASTQs
from collections import Counter

FLANK_LENGTH = 10

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
    min_align_score: float = 0.5,
    mismatches_tolerated: int = 1,
):
    """Extracts TuBa-seq double barcodes from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        FASTQ_directory (or file) containing the forward FASTQ run(s)
        sgRNA_file: All sgRNAs used and their corresponding genes
        master_read: Flanking amplicon sequence to align to each read
        min_align_score: Combined PE alignment score needed to use read, Range [0, 1)
        mismatches_tolerated: # of mismatches tolerated in sgID"""

    # in case running interactively, convert strings to paths
    FASTQ_directory = Path(FASTQ_directory)
    sgRNA_file = Path(sgRNA_file)

    sg_info = pd.read_csv(sgRNA_file, converters={"ID": str.upper})

    sgID_lengths = sg_info["ID"].str.len()
    sgID_length = sgID_lengths[0]
    assert (sgID_lengths == sgID_length).all(), "sgIDs must be equal length."

    sgRNA_map = BarcodeSet(
        sg_info.set_index("ID")["target"],
        n_mismatches=mismatches_tolerated,
        robust=True,
    )

    sgRNA_map.pop_nonunique()

    master_read = MasterRead(master_read)
    max_scores = master_read.max_simplex_scores

    pileups = Counter()
    scores = Counter()

    min_int_score = int(min_align_score * max_scores[-1])

    file_pair = get_PE_FASTQs(FASTQ_directory)
    FASTQ_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in FASTQ_iter:

        simplex_alignment = master_read.simplex_align(fwd_dna, rev_dna)
        score_set = simplex_alignment.get_scores()
        scores[score_set] += 1

        if score_set[-1] <= min_int_score:
            continue

        barcode = simplex_alignment.extract_barcode()

        known_barcode = barcode[:sgID_length]
        sgRNA_target = sgRNA_map.get(known_barcode, "Unknown Target")
        random_barcode = barcode[sgID_length:]

        pileups[(sgRNA_target, random_barcode)] += 1

    ## Output

    scores = pd.Series(
        {
            tuple((idx / max_i for idx, max_i in zip(idxs, max_scores))): counts
            for idxs, counts in scores.items()
        },
        name="reads",
    )

    scores.index.names = ["Fwd-Ref", "Rev-Ref", "Fwd-Rev", "Core-Ref"]

    scores.to_csv(FASTQ_directory / "scores.csv")

    pileups = pd.Series(pileups, name="reads")
    pileups.index.names = "target", "barcode"
    pileups.to_csv(FASTQ_directory / "pileups.csv")

    final_marginal = scores.groupby(level="Core-Ref").sum()
    outcomes = pd.Series(
        {
            "Poor Alignment": final_marginal.loc[:min_align_score].sum(),
            "Index Mismatches": FASTQ_iter.index_mismatches,
            "Passed": pileups.sum(),
            "Unknown Target": pileups.loc["Unknown Target"].sum(),
        },
        name="reads",
    )
    outcomes.index.names = ["outcome"]
    outcomes.to_csv(FASTQ_directory / "outcomes.csv")


if __name__ == "__main__":
    CLI()
