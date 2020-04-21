#! /usr/bin/env python3

import sys

import numpy
from numpy.lib.stride_tricks import as_strided

window_size = 7
min_site_dp = 10
min_out_masked = 10
fold_diff = 5


def window(
    array:numpy.ndarray, window:int,
    shift:int=1, copy:bool=False
):
    shape = (array.size - window + 1, window)
    stride = array.strides * 2
    
    view = as_strided(
        array, strides=stride, shape=shape
    )[0::shift]

    if copy:
        return view.copy()

    else:
        return view


from ace import AceFile, Contig

acefile = AceFile(sys.argv[1])

for x in range(25):
    y = next(acefile)

# beginning preliminary analysis of this contig
# should be the one from VUNXX_CL0260, contig25
unmasked_win_avg = window(y.unmasked, window_size, 7).sum(axis=1) / window_size
masked_win_avg = window(y.masked, window_size, 7).sum(axis=1) / window_size

masked_ratios = masked_win_avg / (unmasked_win_avg + masked_win_avg)
mrs = masked_ratios[1:]
mro = masked_ratios[:-1]

left = (mro >= fold_diff * mrs).argmax() * (window_size + 1)

mro, mrs = mrs[::-1], mro[::-1]
right = -(mro >= fold_diff * mrs).argmax() * (window_size + 1) - 1