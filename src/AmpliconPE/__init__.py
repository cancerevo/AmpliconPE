#!/usr/bin/env python3
from .ssw_lib import SW
import numpy as np
from datetime import datetime


ALIGNMENT_PARAMS = dict(match=3, mismatch=-1, gap_open=6, gap_extend=1)


reverse_map = np.zeros(256, dtype=np.uint8)
for nuc, comp in zip(b"ATCGN", b"TAGCN"):
    reverse_map[nuc] = comp


def reverse_compliment(s):
    return b"".join(reverse_map[np.frombuffer(s, np.uint8)[::-1]])


def smart_open(filename, mode="rb", makedirs=False):
    """Infers compression of file from extension.

    Parameters:
    -----------
    mode : Filemode string (default: 'rb').

    makedirs : Create directory tree for file, if non-existent (default: False)."""
    from pathlib import Path
    import gzip, bz2, lzma

    file_openers = dict(gz=gzip, gzip=gzip, lzma=lzma, xz=lzma, bz2=bz2)
    File = Path(filename)
    if makedirs:
        File.parent.mkdir(exist_ok=True)
    compression = File.suffix[1:]
    if compression in file_openers:
        open = file_openers[compression].open
    return open(str(filename), mode)


def get_PE_FASTQs(directory, read_pattern=r"_R[12]_", fastq_pattern="*.fastq*"):
    import pathlib, re

    fastqs = list(pathlib.Path(directory).glob(fastq_pattern))
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
        if match.group() == "_R1_":
            fwd = i

    return fastqs[fwd], fastqs[1 - fwd]


class pairedFASTQiter(object):
    def __iter__(self):
        return self

    def __init__(self, fwd_file, rev_file, check_indecies=True):
        self.fwd = smart_open(fwd_file)
        self.rev = smart_open(rev_file)
        self.check_indecies = check_indecies
        self.filtered_reads = 0

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

        ILLUMINA_FAILED_FILTER = b":Y:"
        if (
            not ILLUMINA_FAILED_FILTER in fwd_header
            and not ILLUMINA_FAILED_FILTER in rev_header
        ):
            if (
                not self.check_indecies
                or fwd_header.rpartition(b":")[2] == rev_header.rpartition(b":")[2]
            ):
                return fwd_dna, rev_dna
        self.filtered_reads += 1


# See https://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
from itertools import combinations, product


def mismatcher(word, i=2):
    for d in range(i + 1):
        for locs in combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            # try:
            for poss in product(*thisWord):
                yield "".join(poss)
            # except:
            #    for poss in product(*thisWord):
            #        print(poss)
            #    assert False


# See https://realpython.com/inherit-python-dict/
class BarcodeSet(dict):
    def __setitem__(self, barcode, target):
        for mismatch in mismatcher(barcode, self.n_mismatches):
            if mismatch in self:
                if self.robust:
                    super().__setitem__(mismatch, list(self[mismatch]) + [target])
                else:
                    raise ValueError(
                        f"""{mismatch}, a mismatch of {barcode}, is already in BarcodeSet.
There are currently {len(self)} barcodes in this set.
Use BarcodeSet(robust=True) if you would like non-unique mismatches to map to a list of all possible labels."""
                    )

            super().__setitem__(mismatch, target)
        super().__setitem__(barcode, target)

    def update(self, other):
        for k, v in other.items():
            self[k] = v

    def __init__(self, *args, n_mismatches=1, robust=False):
        self.n_mismatches = n_mismatches
        self.robust = robust
        if len(args) == 1:
            self.update(dict(args[0]))


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
        print(fwd_bc)
        print(rev_bc)
        print()
        print(self.fwd_align.nScore, self.rev_align.nScore)
        assert False

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
