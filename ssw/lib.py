
import os
import inspect
import ctypes
import numpy

from kit.utils import rc


bases = ['A', 'C', 'G', 'T', 'N']
base_to_int = {x: i for i, x in enumerate(bases)}
int_to_base = {i: x for i, x in enumerate(bases)}


default_alignment_parameters = {
    'match': 2, 'mismatch': 2, 'gap_open': 3, 'gap_extend': 1, 'nflag': 2
}

cigar_info = 'MIDNSHP=X'


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


class CSmithWaterman:

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

        self.parameters = {
            'match': 0, 'mismatch': 0, 'gap_open': 0, 'gap_extend': 0, 'nflag': 0
        }

        # adding param attrs here, because otherwise editor complains about not having them
        self.match = None
        self.mismatch = None
        self.gap_open = None
        self.gap_extend = None
        self.nflag = None

        self.parameters_set = False

    def set_alignment_params(self, **kwargs):
        for k in kwargs:
            if k in self.parameters and not(self.parameters[k]):
                self.__dict__[k] = kwargs[k]
                self.parameters[k] = 1

        # if a particular parameter is not set - use default
        for p in self.parameters:
            if not(self.parameters[p]):
                self.__dict__[p] = default_alignment_parameters[p]
                self.parameters[p] = 1

        self.scoring = build_score_matrix(
            len(bases), self.match, self.mismatch
        )

        self.parameters_set = True

    def view_alignment_params(self):
        for k in self.parameters:
            print(k, self.__dict__[k])

    def reset_parameters(self, **kwargs):
        # unset parameters and remove them as members of object
        if not(kwargs):
            for k in self.parameters:
                self.parameters[k] = 0
                self.__dict__.pop(k)
        else:
            # simply set parameters in kwargs and leave others alone
            for k in kwargs:
                if k in self.parameters:
                    self.__dict__[k] = kwargs[k]

    def align_sequence_pair(self, query, target):
        if not(self.parameters_set):
            print("Please set parameters before trying to align sequences")
            return None

        results_dict  = {}
        results_dict['Query'] = query
        results_dict['Target'] = target

        query = seq_to_int_array(query)
        target = seq_to_int_array(target)
        # assuming we don't know max SW score size -> using 2
        query_profile = self.ssw_init(
            query, ctypes.c_int32(len(query)), self.scoring, len(bases), 2
        )
        mask_length = max(15, len(query) // 2)

        result = self.ssw_align(
            query_profile, target, ctypes.c_int32(len(target)), self.gap_open,
            self.gap_extend, self.nflag, 0, 0, ctypes.c_int32(mask_length)
        )

        for k, _ in result.contents._fields_:
            results_dict[k] = result.contents.__getattribute__(k)

        # add cigar array to results dict
        results_dict['lCigar'] = [
            result.contents.sCigar[x] for x in range(results_dict['nCigarLen'])
        ]

        build_path(results_dict)

        self.align_destroy(result)

        return results_dict


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


def build_path(results:dict):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = 'Q: '
    sA = '   '
    sR = 'T: '
    nQOff = results['nQryBeg']
    nROff = results['nRefBeg']
    for x in results['lCigar']:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += results['Query'][nQOff : nQOff+n]
            sA += ''.join([
                '|' if results['Query'][nQOff+j] == results['Target'][nROff+j] else '*' for j in range(n)
            ])
            sR += results['Target'][nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += results['Query'][nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += results['Target'][nROff : nROff+n]
            nROff += n

    results['sCigar'] = sCigar
    results['sAlignment'] = "\n".join([sQ, sA, sR])
    # return sCigar, sQ, sA, sR
