#! /usr/bin/env python3

import os
import sys

import numpy
from numpy.lib.stride_tricks import as_strided
from matplotlib import pyplot

window_size = 7
min_site_dp = 10
min_out_masked = 10
fold_diff = 5

from ace import AceFile, Contig, window, SwitchpointFinder

# acefile = AceFile(sys.argv[1])
finder = SwitchpointFinder(sys.argv[1], "")
acefile = finder.acefile

for x in range(acefile.ncontigs):
    y = next(acefile)
    if y.nreads / acefile.nreads > 0.01:
        print(y.name, y.average_rd)
        fn = os.path.splitext(sys.argv[1])[0] + f"_{y.name}.png"
        b, e = finder.find_candidates(y)
        fig = y.generate_figure()
        y.add_candidate_switchpoints_to_fig(fig, (b, e))
        fig.savefig(fn)
        pyplot.close(fig)

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