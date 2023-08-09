#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import SimplexMasterRead as MasterRead, pairedFASTQiter, get_PE_FASTQs
from collections import Counter

def derep(
    FASTQ_directory: Path,
    master_read: str,
    trim_fraction: float = 0.6
):
    """Lossless barcode extraction from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        FASTQ_directory (or file) containing the forward FASTQ run(s)
        master_read: Flanking amplicon sequence to align to each read
        trim_fraction: Fraction of barcode flanks to trim before core alignment."""

    # in case running interactively, convert strings to paths
    FASTQ_directory = Path(FASTQ_directory)

    master_read = MasterRead(master_read, trim_fraction = trim_fraction)

    pileups = Counter()
    scores = Counter()


    file_pair = get_PE_FASTQs(FASTQ_directory)
    FASTQ_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in FASTQ_iter:

        alignment = master_read.align(fwd_dna, rev_dna)
        scores[alignment.get_scores()] += 1
        pileups[(alignemt.final_score, alignment.extract_barcode())] += 1

    ## Output
    max_scores = master_read.self_alignment.get_scores()
    
    score_ixs = np.fromiter(scores.keys(), dtype=np.dtype((float, len(max_scores))))
    normalized_scores = score_ixs.div(max_scores, axis=1)

    score_series = pd.Series(scores.values(), 
        index=pd.MultiIndex.from_arrays(normalized_scores, 
            names=master_read.self_alignment.align_pairs),
                name='reads')

    score_series.to_csv(FASTQ_directory / "scores.csv")

    pileups = pd.Series(pileups.values(), name="reads",
        index = pd.MultiIndex.from_tuples(pileups.keys(), names=['score', 'barcode']))
    pileups.to_csv(FASTQ_directory / "pileups.csv")



        # Think about Integer vs. Float for scores & pileups



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
