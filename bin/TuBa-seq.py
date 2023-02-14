#!/usr/bin/env python3
import pandas as pd
from jsonargparse import CLI

from pathlib import Path
from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, get_PE_FASTQs


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
    FASTQ_directory: Path = Path("FASTQs"),
    sgRNA_file: Path = Path("sgRNA_info.csv"),
    output: Path = Path("output.h5"),
    master_read: str = default_master_read,
    min_align_score: float = 0.75,
    mismatches_tolerated: int = 2,
    parallel: bool = False,
):
    """Extracts TuBa-seq double barcodes from Paired-End (PE) FASTQ reads & dereplicates

    Args:
        forward_reads: Directory (or file) containing the forward FASTQ run(s)
        reverse_reads: Directory (or file) containing the reverse FASTQ run(s)
        sgRNA_file: All sgRNAs used and their corresponding genes
        output: Name of HDF5 output store
        master_read: Flanking amplicon sequence to align to each read
        min_align_score: Combined PE alignment score needed to use read, Range [0, 1)
        mismatches_tolerated: # of mismatches tolerated in sgID
        parallel: Multi-Thread operation?"""
    if parallel:
        from pmap import pmap as map

    directories = [name for name in FASTQ_directory.iterdir() if name.is_dir()]

    sg_info = pd.read_csv(sgRNA_file, converters={"ID": str.upper})
    sgID_length = int(sg_info["ID"].str.len().median())

    sgRNA_map = BarcodeSet(
        sg_info.set_index("ID")["target"], n_mismatches=mismatches_tolerated
    )

    master_read = MasterRead(master_read)

    def derep_barcodes(directory):
        import numpy as np
        from collections import Counter

        pileups = Counter()
        scores = pd.DataFrame(
            np.zeros((master_read.max_score + 1, 2), dtype=np.int64),
            index=pd.Index(
                np.linspace(0, 1, num=master_read.max_score + 1), name="Score"
            ),
        )
        min_int_score = int(min_align_score * master_read.max_score)

        file_pair = get_PE_FASTQs(directory)
        fastq_iter = pairedFASTQiter(*file_pair)
        for fwd_dna, rev_dna in FASTQiter:

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

        pileups = pd.Series(pileups)
        pileups.index.names = "target", "barcode"
        return {
            "pileups": pileups,
            "alignment scores": scores,
            "lost reads": pd.Series(
                {
                    "Poor Alignment": poor_alignment,
                    "Length Mismatch": pileups.sum() - poor_alignment,
                    "Index Mismatch": FASTQiter.index_mismatch,
                }
            ),
        }

    N = len(directories)
    outputs = map(derep_barcodes, directories, N * [master_read], N * [sgRNA_map])

    consolidated_outputs = {
        table_name: pd.concat(
            {
                str(sample_name): output[table_name]
                for sample_name, output in zip(directories, outputs)
                if len(output)[table_name] > 0
            },
            names=["Sample"],
        )
        for table_name in outputs[0].keys()
    }

    with pd.HDFStore(output, "w", complevel=9) as store:
        for name, df in consolidated_outputs.items():
            store.put(name, df)


if __name__ == "__main__":
    CLI()
