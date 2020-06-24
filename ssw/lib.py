
import os
import inspect
import ctypes
import numpy

from kit.utils import rc


bases = ['A', 'C', 'G', 'T', 'N']
base_to_int = {x: i for i, x in enumerate(bases)}
int_to_base = {i: x for i, x in enumerate(bases)}


def seq_to_int_array(seq):
    c_seq_array_dec = len(seq) * ctypes.c_int8
    c_seq_array = c_seq_array_dec()

    for i, e in enumerate(seq):
        if e in base_to_int:
            c_seq_array[i] = base_to_int[e]
        else:
            c_seq_array[i] = bases[-1]

    return c_seq_array


def build_score_matrix(n_elements, match_score, mismatch_score):
    score_matrix = numpy.ones(
        (n_elements, n_elements), dtype=numpy.int8
    ) * -mismatch_score
    numpy.fill_diagonal(score_matrix, match_score)
    matrix_in_c = (ctypes.c_int8 * score_matrix.size)()
    matrix_in_c[:] = score_matrix.flatten()
    return matrix_in_c


class CProfile(ctypes.Structure):
    _fields_ = [
        ('pByte', ctypes.POINTER(ctypes.c_int32)),
        ('pWord', ctypes.POINTER(ctypes.c_int32)),
        ('pRead', ctypes.POINTER(ctypes.c_int8)),
        ('pMat', ctypes.POINTER(ctypes.c_int8)),
        ('nReadLen', ctypes.c_int32),
        ('nN', ctypes.c_int32),
        ('nBias', ctypes.c_uint8)
    ]


class CAlignmentResult(ctypes.Structure):
    _fields_ = [
        ('nScore', ctypes.c_uint16),
        ('nScore2', ctypes.c_uint16),
        ('nRefBeg', ctypes.c_int32),
        ('nRefEnd', ctypes.c_int32),
        ('nQryBeg', ctypes.c_int32),
        ('nQryEnd', ctypes.c_int32),
        ('nRefEnd2', ctypes.c_int32),
        ('sCigar', ctypes.POINTER(ctypes.c_uint32)),
        ('nCigarLen', ctypes.c_int32)
    ] 


class CStripedSmithWaterman:

    def __init__(self, path_to_libssw=None, library_filename="libssw.so", debug=True):
        if not(path_to_libssw):
            modframe = inspect.currentframe()
            filename = inspect.getframeinfo(modframe).filename
            path_to_libssw = f"{os.path.split(filename)[0]}/src"

        self.path = path_to_libssw
        self.full_path = os.path.join(self.path, library_filename)

        if os.path.exists(self.full_path):
            self.ssw = ctypes.cdll.LoadLibrary(self.full_path)
            if debug:
                print(f"{self.full_path} successfully loaded.")
        else:
            print(f"{self.full_path} does not exist. Try compiling library.")
            return

        # ssw initialize
        self.ssw_init = self.ssw.ssw_init
        self.ssw_init.argtypes = [
            ctypes.POINTER(ctypes.c_int8),  # pointer to query
            ctypes.c_int32,                 # length of query
            ctypes.POINTER(ctypes.c_int8),  # pointer to score matrix
            ctypes.c_int32,                 # dimension of matrix
            ctypes.c_int8                   # size of SW score (0, 1 or 2)
        ]
        self.ssw_init.restype = ctypes.POINTER(CProfile)

        # ssw initialize destroy (deallocate)
        self.init_destroy = self.ssw.init_destroy
        self.init_destroy.argtypes = [ctypes.POINTER(CProfile)]
        # memory deallocation returns nothing (void)
        self.init_destroy.restype = None

        # ssw align
        self.ssw_align = self.ssw.ssw_align
        self.ssw_align.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_int8),
            ctypes.c_int32,
            ctypes.c_uint8,
            ctypes.c_uint8,
            ctypes.c_uint8,
            ctypes.c_uint16,
            ctypes.c_int32,
            ctypes.c_int32
        ]
        self.ssw_align.restype = ctypes.POINTER(CAlignmentResult)

        # ssw align destroy
        self.align_destroy = self.ssw.align_destroy
        self.align_destroy.argtypes = [ctypes.POINTER(CAlignmentResult)]
        self.align_destroy.restype = None


def try_alignment(path, query, target):
    match, mismatch = 2, 2
    mat = build_score_matrix(len(bases), match, mismatch)
    ssw = CStripedSmithWaterman(path)
    q = seq_to_int_array(query)
    t = seq_to_int_array(target)

    qprof = ssw.ssw_init(
        q, ctypes.c_int32(len(q)), mat, len(bases), 2
    )
    gap_open = 3
    gap_extend = 1
    nflag = 2
    mask_length = len(query) // 2 if len(query) >= 30 else 15

    res = ssw.ssw_align(
        qprof, t, ctypes.c_int32(len(t)), gap_open, gap_extend, nflag, 0, 0, mask_length
    )

    return res
