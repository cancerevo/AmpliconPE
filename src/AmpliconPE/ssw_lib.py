#!/usr/bin/env python3
"""
Simple python wrapper for SSW library
Please put the path of libssw.so into LD_LIBRARY_PATH or pass it explicitly as a parameter
By Yongan Zhao (March 2016)
"""

import ctypes as ct
import numpy as np

###############################################################################
## load libssw (Chris' extension - very buggy)
###############################################################################


def get_libssw_path():
    from importlib.util import find_spec

    libssw_path = find_spec("libssw")
    if libssw_path is not None:
        return libssw_path.origin
    try:
        from pathlib import Path

        # print(find_spec("AmpliconPE"))
        p = Path(find_spec("AmpliconPE").submodule_search_locations[0])
        return next(p.parent.glob("*libssw*"))
    except:
        from glob import glob

        base_dir = p.parent.parent.parent
        print(base_dir)
        libssw_path = glob(str(base_dir) + "/**/libssw*.so", recursive=True)
        if len(libssw_path):
            return libssw_path[0]
        # return glob("/home/runner/work/AmpliconPE/**/libssw*.so", recursive=True)[0]
        raise ImportError("Could not the libssw shared-object library.")


c_extension = ct.cdll.LoadLibrary(get_libssw_path())

o_N = ord(b"N")
o_gap = ord(b"-")

###############################################################################
###############################################################################


def buffer_merge(qry, ref):
    """buffer_merge(qry, ref) -> Performant consensus sequence of equal-length byte-strings with mismatches replaced by 'N'"""
    if qry == ref:
        return qry
    qry_buffer = np.frombuffer(qry, dtype=np.uint8)
    ref_buffer = np.frombuffer(ref, dtype=np.uint8)
    return np.where(qry_buffer != ref_buffer, o_N, qry_buffer).tobytes()


class CAlignRes(ct.Structure):
    """
    @typedef	structure of the alignment result
    @field	nScore	    the best alignment score
    @field	nScore2	    sub-optimal alignment score
    @field	nRefBeg	    0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning
                                            position is not available
    @field	nRefEnd	    0-based best alignment ending position on reference
    @field	nQryBeg	    0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
                                            position is not available
    @field	nQryEnd	    0-based best alignment ending position on read
    @field	nRefEnd2    0-based sub-optimal alignment ending position on read
    @field	sCigar	    best alignment cigar; stored the same as that in BAM format, high 28 bits: length, low 4 bits: M/I/D (0/1/2);
                                    cigar = 0 when the best alignment path is not available
    @field	nCigarLen   length of the cigar string; cigarLen = 0 when the best alignment path is not available
    """

    _fields_ = [
        ("nScore", ct.c_uint16),
        ("nScore2", ct.c_uint16),
        ("nRefBeg", ct.c_int32),
        ("nRefEnd", ct.c_int32),
        ("nQryBeg", ct.c_int32),
        ("nQryEnd", ct.c_int32),
        ("nRefEnd2", ct.c_int32),
        ("sCigar", ct.POINTER(ct.c_uint32)),
        ("nCigarLen", ct.c_int32),
    ]

    def __dalloc__(self):
        ssw.align_destroy(self)

    def build_cigar(self, query, reference):
        """
        build cigar string and align path based on cigar array returned by ssw_align
        @param  q   query sequence
        @param  r   reference sequence
        @param  nQryBeg   begin position of query sequence
        @param  nRefBeg   begin position of reference sequence
        @param  lCigar   cigar array
        """
        sCigarInfo = b"MIDNSHP=X"
        sQ = []
        sR = []
        nQOff = self.nQryBeg
        nROff = self.nRefBeg
        for idx in range(self.nCigarLen):
            x = self.sCigar[idx]
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

    def build_consensus(self, query, reference, expected_length):
        """
        Build consensus sequence from query/reference pair.

        Replaces mismatch bases with 'N', and includes Insertions/Deletions if they
        bring the consensus sequence closer to the `expected_length`.
        """
        sCigarInfo = b"MIDNSHP=X"
        consensus = []
        nQOff = self.nQryBeg
        nROff = self.nRefBeg

        print("here")

        max_length = sum([self.sCigar[idx] >> 4 for idx in range(self.nCigarLen)])
        too_long = max_length > expected_length

        print(self.sCigar)
        print(max_length)
        print(too_long)

        for idx in range(self.nCigarLen):
            x = self.sCigar[idx]
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
        dBeg = self.nQryBeg - self.nRefBeg
        truncation = min(self.nQryBeg, self.nRefBeg)

        dRef_pre_start = 0

        ref_start_remaining = barcode_start - dBeg - truncation
        idx = 0
        while ref_start_remaining > 0:
            x = self.sCigar[idx]
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
        while ref_barcode_remaining > 0 and idx < self.nCigarLen:
            x = self.sCigar[idx]
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

        return slice(
            barcode_start + dBeg + dRef_pre_start,
            barcode_stop + dBeg + dRef_barcode + dRef_pre_start,
        )

    def query_core(self):
        return slice(self.nQryBeg, self.nQryEnd + 1)


class CProfile(ct.Structure):
    """
    @typedef	structure of the query profile
    @field	pByte	byte array for profile
    @field	pWord	word array for profile
    @field	pRead	number array for read
    @field	pMat	score matrix
    @field	nReadLen	read length
    @field	nN	edge length of score matrix
    @field	nBias	bias
    """

    _fields_ = [
        ("pByte", ct.POINTER(ct.c_int32)),
        ("pWord", ct.POINTER(ct.c_int32)),
        ("pRead", ct.POINTER(ct.c_int8)),
        ("pMat", ct.POINTER(ct.c_int8)),
        ("nReadLen", ct.c_int32),
        ("nN", ct.c_int32),
        ("nBias", ct.c_uint8),
    ]

    def __dalloc__(self):
        ssw.align_destroy(self)


class CSsw(object):
    """
    A class for libssw
    """

    def __init__(self):
        """
        init all para
        @para   sLibpath    argparse object
        @function	Create the query profile using the query sequence.
        @param	read	pointer to the query sequence; the query sequence needs to be numbers
        @param	readLen	length of the query sequence
        @param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
        @param	n	the square root of the number of elements in mat (mat has n*n elements)
        @param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
                                                your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
        @return	pointer to the query profile structure
        @note	example for parameter read and mat:
                        If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
                        Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
                        //A  C  G  T
                          2 -2 -2 -2 //A
                         -2  2 -2 -2 //C
                         -2 -2  2 -2 //G
                         -2 -2 -2  2 //T
                        mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
        """
        self.ssw = c_extension
        self.ssw_init = self.ssw.ssw_init
        self.ssw_init.argtypes = [
            ct.POINTER(ct.c_int8),
            ct.c_int32,
            ct.POINTER(ct.c_int8),
            ct.c_int32,
            ct.c_int8,
        ]
        self.ssw_init.restype = ct.POINTER(CProfile)
        # init init_destroy
        """
	@function	Release the memory allocated by function ssw_init.
	@param	p	pointer to the query profile structure
        """
        self.init_destroy = self.ssw.init_destroy
        self.init_destroy.argtypes = [ct.POINTER(CProfile)]
        self.init_destroy.restype = None
        # init ssw_align
        """
!	@function	Do Striped Smith-Waterman alignment.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
				function ssw_init
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty
	@param	weight_gapE	the absolute value of gap extension penalty
	@param	flag	bitwise FLAG; (from high to low) bit 5: when setted as 1, function ssw_align will return the best alignment
					beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
					< filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
					cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
  					will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
 					setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
					the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
	@param	filters	score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
 					decription of the flag parameter for detailed usage.)
	@param	filterd	distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
					the decription of the flag parameter for detailed usage.)
	@param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
					readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
					return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
					alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
					largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
					picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
					reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
					score from the unmasked elements.
	@return	pointer to the alignment result structure
	@note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
			and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
			while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
			0-based coordinate.
        """
        self.ssw_align = self.ssw.ssw_align
        self.ssw_align.argtypes = [
            ct.c_void_p,
            ct.POINTER(ct.c_int8),
            ct.c_int32,
            ct.c_uint8,
            ct.c_uint8,
            ct.c_uint8,
            ct.c_uint16,
            ct.c_int32,
            ct.c_int32,
        ]
        self.ssw_align.restype = ct.POINTER(CAlignRes)
        # init align_destroy
        """
	@function	Release the memory allocated by function ssw_align.
	@param	a	pointer to the alignment result structure
        """
        self.align_destroy = self.ssw.align_destroy
        self.align_destroy.argtypes = [ct.POINTER(CAlignRes)]
        self.align_destroy.restype = None


ssw = CSsw()


class SW(object):
    def __dealloc__(self):
        ssw.init_destroy(self.qProfile)

    def __del__(self):
        ssw.init_destroy(self.qProfile)

    def __init__(self, match=4, mismatch=2, gap_open=6, gap_extend=1):
        import numpy as np

        # init DNA score matrix
        lEle = list(b"ACGTN")
        nInt2Ele = np.array(lEle, dtype="int8")
        self.nEle2Int = np.zeros(256, dtype="int8")
        self.nEle2Int[nInt2Ele] = np.arange(len(lEle))
        self.lScore = [
            0 if nuc1 == o_N or nuc2 == o_N else (match if nuc1 == nuc2 else mismatch)
            for nuc1 in lEle
            for nuc2 in lEle
        ]
        # translate score matrix to ctypes
        self.mat = (len(self.lScore) * ct.c_int8)()
        self.mat[:] = self.lScore
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.lEle = lEle

    def align(self, query, reference):
        """Performs the underling alignment"""

        nQuery = self.nEle2Int[list(query)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.qProfile = ssw.ssw_init(
            nQuery, ct.c_int32(len(query)), self.mat, len(self.lEle), 2
        )
        nMaskLen = len(query) // 2 if len(query) > 30 else 15
        nFlag = 2
        nReference = self.nEle2Int[list(reference)].ctypes.data_as(
            ct.POINTER(ct.c_int8)
        )
        return ssw.ssw_align(
            self.qProfile,
            nReference,
            ct.c_int32(len(reference)),
            self.gap_open,
            self.gap_extend,
            nFlag,
            0,
            0,
            nMaskLen,
        ).contents
