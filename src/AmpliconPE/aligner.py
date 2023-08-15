#!/usr/bin/env python3
"""
Extension of ssw_lib.py
"""

import ctypes as ct
import numpy as np
from memory_profiler import profile
from .ssw_lib import CSsw

N_memtest = int(1e4)
VERBOSE = False

o_N = ord(b"N")


def buffer_merge(qry, ref):
    """buffer_merge(qry, ref) -> Performant consensus sequence of equal-length byte-strings with mismatches replaced by 'N'"""
    if qry == ref:
        return qry
    qry_buffer = np.frombuffer(qry, dtype=np.uint8)
    ref_buffer = np.frombuffer(ref, dtype=np.uint8)
    return np.where(qry_buffer != ref_buffer, o_N, qry_buffer).tobytes()


def get_libssw_path(verbose=False):
    """get_libssw_path() -> path to libssw.so (Compiled Striped Smith Waterman (SSW) C Library)

    Can be in sys.modules, or in the parent directory containing AmpliconPE (searches for both)."""
    from importlib.util import find_spec
    from pathlib import Path
    from glob import glob

    modules_path = find_spec("libssw")
    parent_dir_search = Path(find_spec("AmpliconPE").submodule_search_locations[0])

    if modules_path is not None:
        directory = "sys.modules"
        path = modules_path.origin
    else:
        try:
            directory = "AmpliconPE parent directory"
            path = next(parent_dir_search)
        except:
            parent_parent_dir = parent_dir_search.parent.parent.parent
            libssw_path = glob(
                str(parent_parent_dir) + "/**/libssw*.so", recursive=True
            )
            try:
                path = libssw_path[0]
                directory = "AmpliconPE parent, parent directory"
            except:
                raise ImportError("Could not find the libssw.so shared-object library.")

    if verbose:
        print(f"libssw.so in {directory} @ {path}")
    return path


c_extension = ct.cdll.LoadLibrary(get_libssw_path(verbose=VERBOSE))
SSW = CSsw(c_extension)


class Alignment(object):
    def __del__(self):
        SSW.align_destroy(self.c_align_res)

    def __init__(self, c_align_res, query, reference):
        self.c_align_res = c_align_res
        self.query = query
        self.reference = reference
        self.score = c_align_res.nScore

    def __str__(self):
        return "{:}\n{:}  Score: {score}".format(*self.build_cigar(), score=self.score)

    def build_cigar(self):
        """
        build cigar string and align path based on cigar array returned by ssw_align
        @param  q   query sequence
        @param  r   reference sequence
        @param  nQryBeg   begin position of query sequence
        @param  nRefBeg   begin position of reference sequence
        @param  lCigar   cigar array
        """
        c_align = self.c_align_res
        sCigarInfo = b"MIDNSHP=X"
        sQ = []
        sR = []
        nQOff = c_align.nQryBeg
        nROff = c_align.nRefBeg

        query = self.query
        reference = self.reference

        for idx in range(c_align.nCigarLen):
            x = c_align.sCigar[idx]
            n = x >> 4
            m = x & 15
            c = 77 if m > 8 else sCigarInfo[m]
            if c == 77:  #'M'
                sQ.append(query[nQOff : nQOff + n])
                sR.append(reference[nROff : nROff + n])
                nQOff += n
                nROff += n
            elif c == 73:  #'I'
                sQ.append(query[nQOff : nQOff + n])
                sR.append(b"-" * n)
                nQOff += n
            elif c == 68:  #'D'
                sQ.append(b"-" * n)
                sR.append(reference[nROff : nROff + n])
                nROff += n
            else:
                raise ValueError("Invalid Cigar Annotation ({:})".format(c))

        return b"".join(sQ), b"".join(sR)

    def build_consensus(self, expected_length):
        """
        Build consensus sequence from query/reference pair.

        Replaces mismatch bases with 'N', and includes Insertions/Deletions if they
        bring the consensus sequence closer to the `expected_length`.
        """
        query = self.query
        reference = self.reference

        sCigarInfo = b"MIDNSHP=X"
        consensus = []

        c_align = self.c_align_res
        nQOff = c_align.nQryBeg
        nROff = c_align.nRefBeg

        max_length = sum([c_align.sCigar[idx] >> 4 for idx in range(c_align.nCigarLen)])
        too_long = max_length > expected_length

        for idx in range(c_align.nCigarLen):
            x = c_align.sCigar[idx]
            n = x >> 4
            m = x & 15
            c = 77 if m > 8 else sCigarInfo[m]

            if c == 77:  #'M' = match
                consensus.append(
                    buffer_merge(query[nQOff : nQOff + n], reference[nROff : nROff + n])
                )
                nQOff += n
                nROff += n
            elif c == 73:  #'I' = query insertion
                if too_long:
                    max_length -= n
                    too_long = max_length > expected_length
                else:
                    consensus.append(query[nQOff : nQOff + n])
                nQOff += n
            elif c == 68:  #'D' = query deletion
                if too_long:
                    max_length -= n
                    too_long = max_length > expected_length
                else:
                    consensus.append(reference[nROff : nROff + n])
                nROff += n
            else:
                raise ValueError("Invalid Cigar Annotation ({:})".format(c))

        return b"".join(consensus)

    def extract_barcode(self, barcode_start, barcode_stop):
        sCigarInfo = b"MIDNSHP=X"
        c_align = self.c_align_res
        dBeg = c_align.nQryBeg - c_align.nRefBeg
        truncation = min(c_align.nQryBeg, c_align.nRefBeg)

        dRef_pre_start = 0

        ref_start_remaining = barcode_start - dBeg - truncation
        idx = 0
        while ref_start_remaining > 0:
            x = c_align.sCigar[idx]
            n = x >> 4
            m = x & 15
            c = 77 if m > 8 else sCigarInfo[m]
            if c == 77:  # "M" = match
                ref_start_remaining -= n
            elif c == 73:  #'I' = Ins
                dRef_pre_start += n
            elif c == 68:  #'D' = del
                ref_start_remaining -= n
            idx += 1

        ref_barcode_remaining = barcode_stop + ref_start_remaining - barcode_start
        dRef_barcode = 0
        while ref_barcode_remaining > 0 and idx < c_align.nCigarLen:
            x = c_align.sCigar[idx]
            n = x >> 4
            m = x & 15
            c = 77 if m > 8 else sCigarInfo[m]
            if c == 77:  # "M" = match
                ref_barcode_remaining -= n
            elif c == 73:  #'I' = Ins
                dRef_barcode += n
            elif c == 68:  #'D' = del
                ref_barcode_remaining -= n
                dRef_barcode -= n
            idx += 1

        return self.query[
            barcode_start
            + dBeg
            + dRef_pre_start : barcode_stop
            + dBeg
            + dRef_barcode
            + dRef_pre_start
        ]

    def query_core(self):
        return self.query[self.c_align_res.nQryBeg : self.c_align_res.nQryEnd + 1]


class Aligner(object):
    alphabet = b"ACGTN"
    gap_open = 6
    gap_extend = 1
    match = 4
    mismatch = 2

    def __init__(self, **kwargs):
        from itertools import product

        self.trans_table = bytes.maketrans(
            self.alphabet, bytes(list(range(len(self.alphabet))))
        )

        self.score_array = [
            0
            if n_1 == o_N or n_2 == o_N
            else (self.match if n_1 == n_2 else self.mismatch)
            for n_1, n_2 in product(self.alphabet, repeat=2)
        ]
        self.mat = (len(self.score_array) * ct.c_int8)()
        self.mat[:] = self.score_array

    def align(self, query, reference):
        """Performs the underling alignment"""

        trans_query = (len(query) * ct.c_int8)()
        trans_query[:] = list(query.translate(self.trans_table))
        trans_reference = (len(reference) * ct.c_int8)()
        trans_reference[:] = list(reference.translate(self.trans_table))

        qProfile = SSW.ssw_init(
            trans_query, ct.c_int32(len(query)), self.mat, len(self.alphabet), 2
        )
        nMaskLen = len(query) // 2 if len(query) > 30 else 15
        nFlag = 2
        alignment = SSW.ssw_align(
            qProfile,
            trans_reference,
            ct.c_int32(len(reference)),
            self.gap_open,
            self.gap_extend,
            nFlag,
            0,
            0,
            nMaskLen,
        )
        SSW.init_destroy(qProfile)
        return Alignment(alignment.contents, query, reference)


@profile
def CSsw_memtest():
    from ssw_lib import CSsw

    ssw = CSsw(c_extension)
    for i in range(N_memtest):
        ssw = CSsw(c_extension)
    return ssw


@profile
def CAlignRes_memtest(read_len=32):
    read = (read_len // 4) * b"ATCG"
    sm = Aligner(gap_open=3)
    read_ins = read[: read_len // 2] + b"C" + read[read_len // 2 :]
    for i in range(N_memtest):
        alignment = sm.align(read_ins, read)
    return alignment.build_cigar()


if __name__ == "__main__":
    CSsw_memtest()
    cigar_tuple = CAlignRes_memtest()
    print(cigar_tuple[0].decode("ascii"))
    print(cigar_tuple[1].decode("ascii"))
