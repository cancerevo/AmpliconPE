#!/usr/bin/env python3
from ssw_lib import SW
import numpy as np
import datetime


ALIGNMENT_PARAMS = dict(
    match=3, 
    mismatch=-1,
    gap_open=6,
    gap_extend=1)


reverse_map = np.zeros(256, dtype=np.uint8)
for nuc, comp in zip(b'ATCGN', b'TAGCN'):
    reverse_map[nuc] = comp

def reverse_compliment(s):
    return b''.join(reverse_map[np.frombuffer(s, np.uint8)[::-1]])
ILLUMINA_FAILED_FILTER = b':Y:'

def smart_open(filename, mode="rb", makedirs=False):
    """Infers compression of file from extension.

    Parameters:
    -----------
    mode : Filemode string (default: 'rb').

    makedirs : Create directory tree for file, if non-existent (default: False)."""
    from pathlib import Path
    import gzip, bz2, lzma
    file_openers = dict(
        gz=gzip, gzip=gzip, lzma=lzma, xz=lzma, bz2=bz2)
    File = Path(filename)
    if makedirs:
        File.parent.mkdir(exist_ok=True)
    compression = File.suffix[1:]
    if compression in file_openers:
        open = file_openers[compression].open 
    return open(str(filename), mode)

class pairedFASTQiter(object):
    def __iter__(self): return self

    def __init__(self, fwd_file, rev_file, check_indecies=True): 
        self.fwd = smart_open(fwd_file)
        self.rev = smart_open(rev_file)
        self.check_indecies = check_indecies
        self.filtered_reads = 0

    def __next__(self):
        fwd_header = self.fwd.readline()
        fwd_dna = self.fwd.readline()
        self.fwd.readline()
        fwd_QC = self.fwd.readline()
        
        rev_header = self.rev.readline()
        rev_dna = self.rev.readline()
        self.rev.readline()
        rev_QC = self.rev.readline()

        if not (fwd_header and rev_header and fwd_QC and rev_QC):
            
            self.fwd.close()
            self.rev.close()
            if not (fwd_header and rev_header):
                raise StopIteration
            if not fwd_QC or not rev_QC:
                raise RuntimeError("Forward & Reverse FASTQ files have different lengths.")
            raise RuntimeError("Input FASTQ file lengths are not a multiple of 4.")

        if not ILLUMINA_FAILED_FILTER in fwd_header and not ILLUMINA_FAILED_FILTER in rev_header:
            if not self.check_indecies or fwd_header.rpartition(b':')[2] == rev_header.rpartition(b':')[2]:
                return fwd_dna, rev_dna
        self.filtered_reads += 1

# See https://stackoverflow.com/questions/11679855/introducing-mutations-in-a-dna-string-in-python
from itertools import combinations,product
def mismatcher(word, i = 2):
    for d in range(i+1):
        for locs in combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in product(*thisWord):
                yield "".join(poss)    

# See https://realpython.com/inherit-python-dict/
class BarcodeSet(dict):
    def __setitem__(self, barcode, target):
        for mismatch in mismatcher(barcode, self.n_mismatches):
            if mismatch in self:
                raise ValueError(f"{mismatch}, a mismatch of {barcode}, is already in BarcodeSet.")
            
            super().__setitem__(mismatch, target)
        super().__setitem__(barcode, target)

    def update(self, other):
        for k, v in other.items(): 
            self[k] = v
    
    def __init__(self, *args, n_mismatches=1):
        self.n_mismatches = n_mismatches
        self.updatae(dict(args))


class DoubleAlignment(object):
    def __init__(self, fwd_align, rev_align, master_read):
        self.fwd_align = fwd_align
        self.rev_align = rev_align
        self.master_read = master_read
        
        self.score = fwd_align.score + rev_align.score

    def extract_barcode(self):
        fwd_bc = self.fwd_align.extract_barcode(self.master_read.barcode_start, self.master_read.barcode_stop)
        rev_bc = reverse_compliment(self.rev_align.extract_barcode(self.master_read.rc_start, self.master_read.rc_stop))

        if len(fwd_bc) != len(rev_bc):
            return "Length Mismatch"
        return fwd_bc if fwd_bc == rev_bc else ''.join([fwd if fwd == rev else 'N' for fwd, rev in zip(fwd_bc, rev_bc)])


class MasterRead(str):
    def align(self, fwd_read, rev_read):
        
        return DoubleAlignment(
            self.sw.align(fwd_read, self.seq),
            self.sw.align(rev_read, self.reverse_compliment),
            self)

    
    def __init__(self, seq, alignment_params=ALIGNMENT_PARAMS):
        self.seq = seq
        self.alignmet_params = alignment_params
        self.barcode_start = seq.index("N")
        self.barcode_stop = seq.rindex("N") + 1
        self.reverse_compliment = reverse_compliment(self.alignment_seq) 
        self.rc_start = self.reverse_compliment.index("N")
        self.rc_stop = self.reverse_compliment.rindex("N") + 1

        self.sw = SW(**self.alignment_params)
        self.self_alignment = self.align(seq, self.reverse_compliment)
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

