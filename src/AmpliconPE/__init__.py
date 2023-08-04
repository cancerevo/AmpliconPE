#!/usr/bin/env python3
from .ssw_lib import SW
import numpy as np
from datetime import datetime


ALIGNMENT_PARAMS = dict(match=3, mismatch=-1, gap_open=6, gap_extend=1)
FASTQ_EXTs = "fastq", "fq"
ILLUMINA_FILTERED = b":Y:"


reverse_map = np.zeros(256, dtype=np.uint8)
for nuc, comp in zip(b"ATCGN", b"TAGCN"):
    reverse_map[nuc] = comp


def reverse_compliment(s):
    return b"".join(reverse_map[np.frombuffer(s, np.uint8)[::-1]])


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
def mismatcher(word, i=2, alterations="ACGTN"):
    """Iterator that yields all possible deviations of `word` up to i differences with `alterations` as all alternate characters."""
    from itertools import combinations, product

    for d in range(i + 1):
        for locs in combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in alterations if l != origChar]
            for poss in product(*thisWord):
                yield "".join(poss)


# See https://realpython.com/inherit-python-dict/
class BarcodeSet(dict):
    def __setitem__(self, barcode, target):
        from scipy.spatial.distance import hamming

        for mismatch in mismatcher(barcode, self.n_mismatches):
            if mismatch in self:
                if self.robust:
                    existing_targets = self[mismatch]
                    super().__setitem__(
                        mismatch,
                        existing_targets | {target}
                        if type(existing_targets) is type(self)
                        else {existing_targets, target},
                    )
                else:
                    other_target = self[mismatch]
                    raise ValueError(
                        f"""
{mismatch}, a {hamming(mismatch, barcode):n}-nt mismatch of {barcode} -> {self[barcode]}, is already in this BarcodeSet, 
as a {hamming(mismatch, self.inverse_base[other_target]):n}-nt mismatch of {self.inverse_base[other_target]} -> {other_target}.
There are currently {len(self)} barcodes in this n_mismatches={self.n_mismatches:n} set.
Use BarcodeSet(robust=True) if you would like non-unique mismatches to map to a set of all possible labels."""
                    )

            super().__setitem__(mismatch, target)
        super().__setitem__(barcode, target)

    def update(self, other):
        for k, v in other.items():
            self.inverse_base[v] = k
            self[k] = v

    def __init__(self, *args, n_mismatches=1, robust=False):
        self.n_mismatches = n_mismatches
        self.robust = robust
        self.inverse_base = dict()
        if len(args) == 1:
            self.update(dict(args[0]))

    def pop_nonunique(self):
        if not self.robust:
            raise RuntimeError(
                "Cannot remove non-unique barcodes from a BarcodeSet that is not `robust`"
            )
        return {k: v for k, v in self.items() if type(v) is type(set)}


class DoubleAlignment(object):
    def __init__(self, fwd_read, rev_read, master_read):
        self.fwd_align = master_read.sw.align(fwd_read, master_read.seq)
        self.rev_align = master_read.sw.align(rev_read, master_read.reverse_compliment)

        self.fwd_read = fwd_read
        self.rev_read = rev_read
        self.master_read = master_read

        self.score = self.fwd_align.nScore + self.rev_align.nScore

    def print_cigars(self):
        fwd_cigar = self.fwd_align.build_cigar(self.fwd_read, self.master_read.seq)
        rev_cigar = self.rev_align.build_cigar(
            self.rev_read, self.master_read.reverse_compliment
        )

        print(
            f"""Fwd Cigar {self.fwd_align.nScore}/{self.master_read.max_score/2}:
--------
{fwd_cigar[0]}
{fwd_cigar[1]}
Reverted Rev Cigar {self.rev_align.nScore}/{self.master_read.max_score/2}:
-----------------------
{reverse_compliment(rev_cigar[0])}
{reverse_compliment(rev_cigar[1])}
"""
        )

    def extract_barcode(self):
        fwd_bc = self.fwd_read[
            self.fwd_align.extract_barcode(
                self.master_read.barcode_start, self.master_read.barcode_stop
            )
        ]
        rev_bc = reverse_compliment(
            self.rev_read[
                self.rev_align.extract_barcode(
                    self.master_read.rc_start, self.master_read.rc_stop
                )
            ]
        )
        if len(fwd_bc) != len(rev_bc):
            return "Length Mismatch"
        if fwd_bc == rev_bc:
            return fwd_bc.decode("ascii")
        fwd_buffer = np.frombuffer(fwd_bc, dtype=np.uint8)
        rev_buffer = np.frombuffer(rev_bc, dtype=np.uint8)

        return (
            np.where(fwd_buffer != rev_buffer, 78, fwd_buffer).tobytes().decode("ascii")
        )  # Performant reassignment of barcode differences to 'N'

    def extract_flanks(self):
        pass

    def flank_mismatches_indels(self):
        pass


class MasterRead(str):
    def align(self, fwd_read, rev_read):
        return DoubleAlignment(fwd_read, rev_read, self)

    def __init__(self, seq, alignment_params=ALIGNMENT_PARAMS):
        self.seq = seq.encode("ascii") if type(seq) == str else seq
        self.alignment_params = alignment_params
        self.barcode_start = self.seq.index(b"N")
        self.barcode_stop = self.seq.rindex(b"N") + 1
        self.reverse_compliment = reverse_compliment(self.seq)
        self.rc_start = self.reverse_compliment.index(b"N")
        self.rc_stop = self.reverse_compliment.rindex(b"N") + 1

        self.sw = SW(**self.alignment_params)
        self.self_alignment = self.align(self.seq, self.reverse_compliment)
        self.max_score = self.self_alignment.score

    def _count_scores(self, double_alignment):
        """Probably deprecate -> easier to count 'N' bases in pileups than to look at a wider fwd-rev match"""
        if not hasattr(self, "scores"):
            self.scores = np.zeros((self.max_score + 1, 3), dtype=np.uint64)

        self.scores[double_alignment.fwd_align.nScore, 0] += 1
        self.scores[double_alignment.rev_align.nScore, 1] += 1
        self_score = self.sw.align(
            double_alignment.fwd_read, reverse_compliment(double_alignment.rev_read)
        ).nScore
        self.scores[self_score, 2] += 1


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
