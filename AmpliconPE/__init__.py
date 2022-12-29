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


class MasterRead(str):
   
    def score(self, FWD_read, REV_read):
        return self.sw.score(self.seq, FWD_read) + self.sw.score(self.reverse_compliment, REV_read)
     
    def __init__(self, seq, alignment_params=ALIGNMENT_PARAMS):
        self.seq = seq
        self.alignmet_params = alignment_params
        self.barcode_start = seq.index("N")
        self.barcode_stop = seq.rindex("N") + 1
        self.reverse_compliment = reverse_compliment(self.alignment_seq) 
        self.rc_start = self.reverse_compliment.index("N")
        self.rc_stop = self.reverse_compliment.rindex("N") + 1

        self.sw = SW(**self.alignment_params)
        self.max_score = self.score(self.alignment_seq, self.reverse_compliment)

    def extract_barcode(self, FWD_read, REV_read):
        fwd_bc = self.sw.extract_barcode(FWD_read, self.seq, self.barcode_start, self.barcode_stop)
        rev_bc = reverse_compliment(self.sw.extract_barcode(REV_read, self.reverse_compliment, self.rc_start, self.rc_stop))
        if len(fwd_bc) != len(rev_bc):
            return "Length Mismatch"
        if fwd_bc == rev_bc:
            return fwd_bc
        return ''.join([fwd if fwd == rev else 'N' for fwd, rev in zip(fwd_bc, rev_bc)])


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

