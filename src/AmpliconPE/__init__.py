#!/usr/bin/env python3
from .aligner import Aligner, buffer_merge
import numpy as np
from datetime import datetime
import pandas as pd

FASTQ_EXTs = "fastq", "fq"
ILLUMINA_FILTERED = b":Y:"

ALIGNMENT_PARAMS = dict(match=2, mismatch=-2, gap_open=6, gap_extend=1)

nucleotides = "ACGTN"


reverse_map = np.zeros(256, dtype=np.uint8)
for nuc, comp in zip(b"ATCGN", b"TAGCN"):
    reverse_map[nuc] = comp


def reverse_compliment(s):
    return b"".join(reverse_map[np.frombuffer(s, np.uint8)[::-1]])


def barcode_content(barcodes):
    return pd.concat(
        {nuc: barcodes.str.count(nuc) for nuc in nucleotides}, names=["content"]
    )


def open_FASTQ(filename):
    """opens potentially-compressed FASTQ properly."""
    from pathlib import Path
    import gzip, bz2, lzma

    File = Path(filename)
    suffix = File.suffix[1:]
    if suffix in FASTQ_EXTs:
        return open(File)

    file_openers = dict(gz=gzip, gzip=gzip, lzma=lzma, xz=lzma, bz2=bz2)
    if suffix not in file_openers:
        raise ValueError(f"Cannot determine compression of {File}.")
    return file_openers[suffix].open(str(filename))


def get_PE_FASTQs(directory, read_pattern=r"_R[12]_"):
    import pathlib, re

    fastqs = [
        filename
        for fastq_ext in FASTQ_EXTs
        for filename in pathlib.Path(directory).glob("*." + fastq_ext + "*")
    ]
    if len(fastqs) != 2:
        raise RuntimeError(
            f"Found {len(fastqs)} fastq files in {directory} (expected 2)"
        )

    for i, fastq in enumerate(fastqs):
        match = re.search(read_pattern, str(fastq))
        if not match.group():
            raise ValueError(
                f'No FASTQ number found in {fastq}, using "{read_pattern}" for pattern search.'
            )
        if "1" in match.group():
            fwd = i

    return fastqs[fwd], fastqs[1 - fwd]


class pairedFASTQiter(object):
    def __iter__(self):
        return self

    def __init__(self, fwd_file, rev_file, check_indecies=True):
        self.fwd = open_FASTQ(fwd_file)
        self.rev = open_FASTQ(rev_file)
        self.check_indecies = check_indecies
        self.index_mismatches = 0

    def __next__(self):
        fwd_header = self.fwd.readline()
        fwd_dna = self.fwd.readline()[:-1]
        self.fwd.readline()
        fwd_QC = self.fwd.readline()[:-1]

        rev_header = self.rev.readline()
        rev_dna = self.rev.readline()[:-1]
        self.rev.readline()
        rev_QC = self.rev.readline()[:-1]

        if not (fwd_header and rev_header and fwd_QC and rev_QC):

            self.fwd.close()
            self.rev.close()
            if not (fwd_header and rev_header):
                raise StopIteration
            if not fwd_QC or not rev_QC:
                raise RuntimeError(
                    "Forward & Reverse FASTQ files have different lengths."
                )
            raise RuntimeError("Input FASTQ file lengths are not a multiple of 4.")

        if not (ILLUMINA_FILTERED in fwd_header or ILLUMINA_FILTERED in rev_header):
            if (
                self.check_indecies
                and fwd_header.rpartition(b":")[2] != rev_header.rpartition(b":")[2]
            ):
                self.index_mismatches += 1
            else:
                return fwd_dna, rev_dna


# See https://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
def mismatcher(word, mismatches, alterations="ACGTN", InDels=True):
    """Iterator that yields all possible deviations of `word` up to i differences with `alterations` as all alternate characters."""
    from itertools import combinations, product

    if mismatches == 0:
        raise ValueError(f"Mismatches must be > 0.")
    if mismatches > len(word):
        raise ValueError("{mismatches=} must be < {len(word)=}.")
    deletions = set()
    for d in range(1, mismatches + 1):
        for locs in combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in alterations if l != origChar]
            for poss in product(*thisWord):
                yield "".join(poss)
            if InDels:  # Deletions
                deletion = "".join((nuc for i, nuc in enumerate(word) if i not in locs))
                if deletion not in deletions:
                    deletions.add(deletion)
                    yield deletion

        if InDels:  # Insertions
            insertions = set()
            for locs in combinations(range(len(word) + 1), d):
                start = word[: locs[0]]
                for inserts in product(alterations, repeat=d):
                    insertion = start + "".join(
                        (insert + word[loc:] for insert, loc in zip(inserts, locs))
                    )
                    if insertion not in insertions:
                        insertions.add(insertion)
                        yield insertion


# See https://realpython.com/inherit-python-dict/
class BarcodeSet(dict):
    def __setitem__(self, barcode, target):
        from scipy.spatial.distance import hamming

        super().__setitem__(barcode, target)
        for mismatch in mismatcher(barcode, self.n_mismatches, InDels=self.InDels):
            if mismatch in self and self[mismatch] != target:
                other_targets = self[mismatch]
                if not self.robust:
                    raise ValueError(
                        f"""
{mismatch}, a {hamming(mismatch, barcode):n}-nt mismatch of {barcode} -> {self[barcode]}, is already in this BarcodeSet, 
as a {hamming(mismatch, self.inverse_base[other_targets]):n}-nt mismatch of {self.inverse_base[other_targets]} -> {other_targets}.
There are currently {len(self)} barcodes in this n_mismatches={self.n_mismatches:n} set.
Use RobustBarcodeSet, if you would like non-unique mismatches to map to a set of all possible labels."""
                    )
                if type(other_targets) is not set:
                    super().__setitem__(mismatch, {other_targets, target})
                else:
                    other_targets |= {target}
                    super().__setitem__(mismatch, other_targets)
            else:
                super().__setitem__(mismatch, target)

    def update(self, other, errors=True):
        for k, v in other.items():
            self.inverse_base[v] = k
            if errors:
                self[k] = v
            else:
                super().__setitem__(k, v)

    def __init__(self, *args, n_mismatches=1, robust=True, InDels=True):
        self.n_mismatches = n_mismatches
        self.InDels = InDels
        self.robust = robust
        self.inverse_base = dict()
        if len(args) == 1:
            self.update(dict(args[0]))

    def pop_nonunique(self):
        if not self.robust:
            raise ValueError("Must be Robust BarcodeSet to have non-unique targets.")
        nonunique = {k for k, v in self.items() if type(v) is set}
        return {k: self.pop(k) for k in nonunique}


class DoubleAlignment(object):
    score_pairs = ["fwd-ref", "rev-ref"]

    def get_scores(self):
        return self.fwd_align.score, self.rev_align.score

    def __init__(self, fwd_read, rev_read, master_read):
        self.fwd_align = master_read.aligner.align(fwd_read, master_read.seq)
        self.rev_align = master_read.aligner.align(
            rev_read, master_read.reverse_compliment
        )

        self.fwd_read = fwd_read
        self.rev_read = rev_read
        self.master_read = master_read
        self.final_score = self.fwd_align.score + self.rev_align.score

    def __str__(self):
        fwd_cigar = self.fwd_align.build_cigar()
        rev_cigar = self.rev_align.build_cigar()

        return f"""Fwd Cigar {self.fwd_align.score}/{self.master_read.max_score/2}:
--------
{fwd_cigar[0]}
{fwd_cigar[1]}
Reverted Rev Cigar {self.rev_align.score}/{self.master_read.max_score/2}:
-----------------------
{reverse_compliment(rev_cigar[0])}
{reverse_compliment(rev_cigar[1])}
"""

    def extract_barcode(self):
        fwd_bc = self.fwd_align.extract_barcode(
            self.master_read.barcode_start, self.master_read.barcode_stop
        )
        rev_bc = reverse_compliment(
            self.rev_align.extract_barcode(
                self.master_read.rc_start, self.master_read.rc_stop
            )
        )
        return (
            buffer_merge(fwd_bc, rev_bc).decode("ascii")
            if len(fwd_bc) == len(rev_bc)
            else "Length Mismatch"
        )


class SimplexAlignment(object):
    align_pairs = ["fwd-ref", "rev-ref", "fwd-rev", "core-ref"]

    def __init__(self, fwd_read, rev_read, master_read):
        self.fwd_align = master_read.aligner.align(fwd_read, master_read.seq)
        self.rev_align = master_read.aligner.align(
            rev_read, master_read.reverse_compliment
        )

        self.fwd_read = fwd_read
        self.rev_read = rev_read
        self.master_read = master_read

        self.fwd_core = self.fwd_align.query_core()
        self.rev_core = reverse_compliment(self.rev_align.query_core())
        self.core_align = master_read.aligner.align(self.fwd_core, self.rev_core)
        self.core_consensus = self.core_align.build_consensus(len(master_read.seq))
        self.final_align = master_read.aligner.align(
            self.core_consensus, master_read.core_seq
        )
        self.final_score = self.final_align.score

    def get_scores(self):
        return (
            self.fwd_align.score,
            self.rev_align.score,
            self.core_align.score,
            self.final_align.score,
        )

    def extract_barcode(self):
        return self.final_align.extract_barcode(
            self.master_read.core_barcode_start, self.master_read.core_barcode_stop
        ).decode("ascii")


class MasterRead(object):
    alignment_params = ALIGNMENT_PARAMS.copy()

    def align(self, fwd_read, rev_read):
        return DoubleAlignment(fwd_read, rev_read, self)

    def __init__(self, seq, trim_fraction=0.75, **alignment_params):
        seq = seq.encode("ascii") if type(seq) == str else seq
        self.seq = seq

        self.alignment_params.update(alignment_params)

        if not b"N" in seq:
            raise ValueError(f"No barcode found in {seq}.")
        self.barcode_start = seq.index(b"N")
        self.barcode_stop = seq.rindex(b"N") + 1
        self.reverse_compliment = reverse_compliment(seq)
        self.rc_start = self.reverse_compliment.index(b"N")
        self.rc_stop = self.reverse_compliment.rindex(b"N") + 1

        self.core_start = int(trim_fraction * self.barcode_start)
        self.core_stop = (
            len(seq) + 1 - int((len(seq) + 1 - self.barcode_stop) * trim_fraction)
        )
        self.core_seq = seq[self.core_start : self.core_stop]

        self.core_barcode_start = self.barcode_start - self.core_start
        self.core_barcode_stop = self.barcode_stop - self.core_start

        self.aligner = Aligner(**self.alignment_params)

        perfect_barcode = self.seq.replace(b"N", b"A")
        self.self_alignment = self.align(
            perfect_barcode, reverse_compliment(perfect_barcode)
        )
        self.max_score = self.self_alignment.final_score


class SimplexMasterRead(MasterRead):
    def align(self, fwd_read, rev_read):
        return SimplexAlignment(fwd_read, rev_read, self)


class logPrint(object):
    def line_break(self):
        self.f.write(80 * "-" + "\n")

    def __call__(self, line, print_line=False, header=False):
        if self.verbose or print_line:
            print(line)
        if header:
            self.f.write((len(line) + 4) * "#" + "\n")
            self.f.write("# " + str(line) + " #\n")
            self.f.write((len(line) + 4) * "#" + "\n")
        else:
            self.f.write(str(line) + "\n")
        self.f.flush()

    def close_logPrint(self):
        runtime = datetime.now() - self.start_time
        self("Runtime: {:}".format(str(runtime).split(".")[0]))
        self.line_break()
        self.f.close()

    def __init__(self, input_args, filename=None):
        import __main__ as main
        import os, atexit

        self.start_time = datetime.now()
        self.program = os.path.basename(main.__file__).partition(".py")[0]
        self.filename = self.program + ".LOG" if filename is None else filename
        args_dict = input_args.__dict__.copy()
        self.verbose = args_dict.pop("verbose", False)
        print("Logging output to", self.filename)
        self.f = open(self.filename, "a")
        self.f.write("\n")
        self.line_break()
        self.f.write(
            "Output Summary of {0.program}, executed at {0.start_time:%c} with the following input arguments:\n".format(
                self
            )
        )
        self.line_break()
        for arg, val in args_dict.items():
            self.f.write("{:}: {:}\n".format(arg, val))
        self.line_break()
        atexit.register(self.close_logPrint)
