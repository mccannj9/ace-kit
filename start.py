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
        fn = os.path.splitext(sys.argv[1])[0] + f"_{y.name}.png"
        candidates, derivatives = finder.find_new_candidates(y)
        fig = y.generate_figure()
        max_depth = (y.unmasked + y.masked).max()

        for i, x in enumerate(candidates):

            if x:
                dx = derivatives[i]
                pos = i - y.shift
                if dx > 0:
                    seq = y.seq[pos:pos+30].replace("*", "-")
                else:
                    seq = y.seq[pos+1-30:pos+1].replace("*", "-")

                print(f">{y.name}_{i}_{pos}_{round(dx)}\n{seq}")

                for j, ax in enumerate(fig.axes):
                    ax.vlines(y.min + i, 0, max_depth, linestyles='dotted')
        fig.savefig(fn)
        pyplot.close(fig)