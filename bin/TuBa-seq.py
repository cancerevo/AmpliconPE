#!/usr/bin/env python3
import pandas as pd
import argparse
from pathlib import Path
from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, logPrint


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

parser = argparse.ArgumentParser(
    description="""Determines tumor number from Paired-End reads.""",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

############################## Input #########################################
IO_group = parser.add_argument_group("IO", "Input/Output Optional Arguments")

IO_group.add_argument(
    "--forward_reads",
    type=Path,
    default=Path("forward_reads"),
    help="Directory (or file) containing the forward FASTQ files (or reads).",
)

IO_group.add_argument(
    "--reverse_reads",
    type=Path,
    default=Path("reverse_reads"),
    help="Directory (or file) containing the reverse FASTQ files (or reads).",
)

IO_group.add_argument(
    "--sgRNA_file",
    type=Path,
    default=Path.home() / "tuba" / "sgRNA_info.csv",
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
    help="Outline of the amplicon sequence, degenerate bases can be either 'N' or '.'. --trim and --symmetric_flanks depend on this being the full length FASTQ sequence that you expect after merging reads.",
)

OP_group = parser.add_argument_group("OP", "Optional arguments affecting operation")

OP_group.add_argument(
    "--min_align_score",
    type=float,
    default=0.75,
    help="Minimum alignment score needed to keep read, Range [0, 1).",
)

OP_group.add_argument(
    "-mismatches_tolerated",
    type=int,
    default=2,
    help="# of mismatches tolerated in sgID",
)

OP_group.add_argument(
    "-p", "--parallel", action="store_true", help="Parallelize operation."
)

args = parser.parse_args()

Log = logPrint(args)
if args.parallel:
    from pmap import pmap as map


if args.forward_reads.is_dir():
    forward_files = list(args.forward_reads.glob("*.fastq*"))
    reverse_files = list(args.reverse_reads.glob("*.fastq*"))
    reverse_basenames = {f.name for f in reverse_files}
    for f in forward_files:
        if f.name not in reverse_basenames:
            raise ValueError(
                "FASTQ filenames in forward and reverse directories do not match."
            )

else:
    forward_files = [args.forward_reads]
    reverse_files = [args.reverse_reads]

sg_info = pd.read_csv(args.sgRNA_file, converters={"ID": str.upper})
sgID_length = int(sg_info["ID"].str.len().median())

sgID_map = BarcodeSet(
    sg_info.set_index("ID")["target"], n_mismatches=args.mismatches_tolerated
)

master_read = MasterRead(args.master_read)


def derep_barcodes(FASTQiter):
    import numpy as np
    from collections import Counter

    pileups = Counter()
    scores = pd.DataFrame(
        np.zeros((master_read.max_score + 1, 2), dtype=np.int64),
        index=pd.Index(np.linspace(0, 1, num=master_read.max_score + 1), name="Score"),
    )
    min_int_score = int(args.min_align_score * master_read.max_score)

    for fwd_dna, rev_dna in FASTQiter:

<<<<<<< HEAD
        double_alignment = master_read.align(fwd_dna, rev_dna)
        scores[double_alignment.score] += 1
        if double_alignment.score <= min_int_score:
            continue  

        barcode = double_alignment.extract_barcode()
        if barcode == 'Length Mismatch':
=======
        score = master_read.score(fwd_dna, rev_dna)
        scores[score] += 1
        if score <= min_int_score:
            continue

        barcode = master_read.extract_barcode(fwd_dna, rev_dna)
        if barcode == "Length Mismatch":
>>>>>>> 8d7446ce195d72789121e36bc06b1085428760c7
            continue

        known_barcode = barcode[:sgID_length]
        sgRNA_target = sgID_map.get(known_barcode, "Unknown Target")
        random_barcode = barcode[sgID_length:]

        pileups[(sgRNA_target, random_barcode)] += 1

    poor_alignment = scores[:min_int_score].sum()
    lost = {
<<<<<<< HEAD
        'Poor Alignment' : poor_alignment,
        'Length Mismatch': pileups.sum() - poor_alignment,
        'Index Mismatch' : FASTQiter.index_mismatch
           }
=======
        "Poor Alignment": poor_alignment,
        "Length Mismatch": pileups.sum() - poor_alignment,
        "Index Mismatch": FASTQiter.index_mismatch,
    }
>>>>>>> 8d7446ce195d72789121e36bc06b1085428760c7
    return pileups, lost, scores


outputs = map(
    derep_barcodes,
    [pairedFASTQiter(*file_pair) for file_pair in zip(forward_files, reverse_files)],
)
sample_names = [f.name.partition(".")[0] for f in forward_files]


output_dfs = [
    pd.concat(
        {
            name: output_S
            for name, output_S in zip(sample_names, output_set)
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
