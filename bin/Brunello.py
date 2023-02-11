#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path
from AmpliconPE import (
    MasterRead,
    BarcodeSet,
    get_PE_FASTQs,
    pairedFASTQiter,
    logPrint,
)


#                AGCTTGTGGAAAGGACGAAACACCG GGGGACGGATCCTGGACATG
#                                                       TTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTATTATAAATCTAAGTCTTTAAA
full_amplicon = (
    "AGCTTGTGGAAAGGACGAAACACCG"
    + 20 * "N"
    + "TTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTATTATAAATCTAAGTCTTTAAA"
)

FLANK_LENGTH = 12

default_master_read = full_amplicon[
    full_amplicon.find("N")
    - FLANK_LENGTH : full_amplicon.rindex("N")
    + 1
    + FLANK_LENGTH
]

parser = argparse.ArgumentParser(
    description="""Determines sgRNA counts from Paired-End Brunello library reads.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

############################## Input #########################################
IO_group = parser.add_argument_group("IO", "Input/Output Optional Arguments")

IO_group.add_argument(
    "--FASTQ_directory",
    type=Path,
    default=Path("FASTQs"),
    help="Iterates over all FASTQs in directory",
)

IO_group.add_argument(
    "--sgRNA_file",
    type=Path,
    default="broadgpp-brunello-library-contents.txt",
    help="All sgRNAs used and their corresponding identifiers.",
)

IO_group.add_argument("-v", "--verbose", action="store_true", help="Output more Info")

IO_group.add_argument(
    "--output", type=Path, default=Path("output.h5"), help="Name of HDF5 output store."
)

IO_group.add_argument(
    "--master_read",
    type=str,
    default=default_master_read,
    help="Non-degenerate amplicon sequence expected",
)

OP_group = parser.add_argument_group("OP", "Optional arguments affecting operation")

OP_group.add_argument(
    "--min_align_score",
    type=float,
    default=0.75,
    help="Minimum alignment score needed to keep read, Range [0, 1).",
)

OP_group.add_argument(
    "--mismatches_tolerated",
    type=int,
    default=1,
    help="# of mismatches tolerated in sgRNA",
)

OP_group.add_argument(
    "-p", "--parallel", action="store_true", help="Parallelize operation."
)

args = parser.parse_args()

Log = logPrint(args)
if args.parallel:
    from pmap import pmap as map


directories = [name for name in args.FASTQ_directory.iterdir() if name.is_dir()]


## Setup sgRNA_map

# Data Table can be found at: https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt
#
# Sample data:
#
# Target Gene ID                                                   1.0                             1.0  ...                    NaN                    NaN
# Target Gene Symbol                                              A1BG                            A1BG  ...  Non-Targeting Control  Non-Targeting Control
# Target Transcript                                        NM_130786.3                     NM_130786.3  ...                    NaN                    NaN
# Genomic Sequence                                        NC_000019.10                    NC_000019.10  ...                    NaN                    NaN
# Position of Base After Cut (1-based)                      58351502.0                      58350637.0  ...                    NaN                    NaN
# Strand                                                         sense                       antisense  ...                    NaN                    NaN
# sgRNA Target Sequence                           CATCTTCTTTCACCTGAACG            CTCCGGGGAGAACTCCGGCG  ...   TTTTTAATACAAGGTAATCT   TTTTTCTCACCCGATGAATC
# Target Context Sequence               ATCGCATCTTCTTTCACCTGAACGCGGTGG  CCGGCTCCGGGGAGAACTCCGGCGCGGGCA  ...                    NaN                    NaN
# PAM Sequence                                                     CGG                             CGG  ...                    NaN                    NaN
# Exon Number                                                      5.0                             6.0  ...                    NaN                    NaN
# Rule Set 2 score                                               0.617

ko_library = pd.read_csv(
    args.sgRNA_file, sep="\t", index_col=["sgRNA Target Sequence"]
)["Target Gene Symbol"]


sgRNA_map = BarcodeSet(
    ko_library.index
    + "_"
    + ko_library,  # Hard to mint simple names...Gene Symbol are non-unique & target sequences are uninformative
    n_mismatches=args.mismatches_tolerated,
    robust=True,
)

overlapping = sum(map(lambda v: len(list(v)) > 1, sgRNA_map.values()))

Log(
    f"{overlapping} of {len(sgRNA_map) - len(ko_library)} sgRNA mismatches overlap ({overlapping/(len(sgRNA_map) - len(ko_library)):.2%})."
)

master_read = MasterRead(args.master_read)


def derep_barcodes(directory):
    import numpy as np
    from collections import Counter

    pileups = Counter()
    scores = np.zeros(master_read.max_score + 1, dtype=np.int64)
    min_int_score = int(args.min_align_score * master_read.max_score)

    file_pair = get_PE_FASTQs(directory)
    fastq_iter = pairedFASTQiter(*file_pair)
    for fwd_dna, rev_dna in fastq_iter:

        double_alignment = master_read.align(fwd_dna, rev_dna)
        double_alignment.print_cigars()
        scores[double_alignment.score] += 1
        if double_alignment.score <= min_int_score:
            continue

        sgRNA = double_alignment.extract_barcode()
        if sgRNA == "Length Mismatch":
            continue

        print(fwd_dna)
        print(rev_dna)
        print(sgRNA)
        # print(f"extracted {sgRNA} -> {sgRNA_map.get(sgRNA, 'Unknown')}")
        pileups[sgRNA_map.get(sgRNA, "Unknown Target")] += 1

    poor_alignment = scores[:min_int_score].sum()
    lost = {
        "Poor Alignment": poor_alignment,
        "Length Mismatch": pileups.sum() - poor_alignment,
        "Index Mismatch": fastq_iter.index_mismatch,
    }
    return pileups, lost, pd.Series(scores)


outputs = map(derep_barcodes, directories)


output_dfs = [
    pd.concat(
        {
            str(name): output_S
            for name, output_S in zip(directories, output_set)
            if len(output_S) > 0
        },
        names=["Sample"],
        axis=0,
    )
    for output_set in zip(*list(outputs))
]


###############################################################################
# Output
###############################################################################

store = pd.HDFStore(args.output, "w", complevel=9)
for name, df in zip(["pileups", "lost", "scores"], output_dfs):
    if name == "pileups":
        df.index.names = "Sample", "target", "barcode"
    store.put(name, df)
store.close()
