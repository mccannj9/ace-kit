
from collections import namedtuple

import numpy

from matplotlib import pyplot

from kit.ace import AceFile
from kit.contig import Contig
from kit.utils import window

Result = namedtuple('Result', ('contig', 'candidates', 'derivatives'))

class SwitchpointFinder:
    def __init__(
        self, input_fn, outdir="./", window_size=7, min_depth=10, min_read_prop=0.1
    ):

        self.acefile = AceFile(input_fn)
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.outdir = outdir
        # self.fasta = open(f"{outdir}/results.fas", "w")
    
    def fit(self):
        contig_dict = {}

        with open(f"{self.outdir}/results.fas", "w") as self.fasta:

            for _ in range(self.acefile.ncontigs):
                ctg = next(self.acefile)

                if ctg.nreads / self.acefile.nreads > self.min_read_prop:
                    cands, derivs = self.find_candidates(ctg)
                    contig_dict[ctg.name] = Result(ctg, cands, derivs)
                    self.write_and_plot_results(contig_dict[ctg.name])

        return contig_dict
    
    def write_and_plot_results(self, result:Result):
        contig = result.contig
        fig = contig.generate_figure()
        max_dp = contig.depth.max()

        for i in numpy.flatnonzero(result.candidates):
            dx = result.derivatives[i]
            pos = i - contig.shift

            if dx > 0:
                seq = contig.seq[pos:pos+30].replace("*", "-")
            else:
                seq = contig.seq[pos+1-30:pos+1].replace("*", "-")
            
            print(f">{contig.name}_{i}_{pos}_{round(dx)}\n{seq}", file=self.fasta)

            for _, ax in enumerate(fig.axes):
                ax.vlines(contig.min + i, 0, max_dp, linestyles='dotted')

        fig.savefig(f"{self.outdir}/{contig.name}.png")
        pyplot.close(fig)

    def find_candidates(self, contig):
        dmask = contig.depth > self.min_depth
        diff = contig.unmasked - contig.masked
        s = numpy.sign(diff)
        sc = ((numpy.roll(s, 1) - s) != 0).astype(int)
        side = self.window_size // 2
        derivatives = numpy.zeros(shape=diff.shape, dtype=float)
        for c in numpy.flatnonzero(sc):
            if (c-side < 0) or (c + side + 1 > contig.depth.size):
                sc[c] = 0
            else:
                wdiff = contig.unmasked[c-side:c+side+1] - contig.masked[c-side:c+side+1]
                wdepth = contig.depth[c-side:c+side+1].mean()
                derivative = (wdiff[-1] - wdiff[0]) / self.window_size
                derivatives[c] = derivative
                if abs(derivative) / wdepth < 0.10:
                    sc[c] = 0
        
        # get the idx of min and max of derivative
        min_der_idx = derivatives.argmin()
        max_der_idx = derivatives.argmax()
        derivatives_mask = numpy.zeros(shape=derivatives.shape)
        derivatives_mask[min_der_idx] = 1
        derivatives_mask[max_der_idx] = 1

        return sc * dmask * derivatives_mask, derivatives