
from collections import namedtuple
from dataclasses import dataclass, field

import numpy

from matplotlib import pyplot

from kit.ace import AceFile
from kit.utils import create_logo, extract_seqs, get_reads_from_candidate
from kit.contig import Contig, Boundary

SingleResult = namedtuple(
    'Result', ('contig', 'position', 'dvalue', 'depth', 'reads')
)


@dataclass
class Result:
    contig: Contig
    candidates: numpy.ndarray
    derivatives: numpy.ndarray
    n: int = 0
    blast_results: list = field(default_factory=list)
    logos: list = field(default_factory=list)

    def __lt__(self, other):
        return self.n < other.n

    def max_mag_derivative(self):
        return numpy.abs(self.derivatives.max())

    def write_contig_boundaries_as_fasta(self, filename:str):
        with open(filename, 'w') as fasta:
            for c in numpy.flatnonzero(self.candidates):
                pos = c - self.contig.shift
                dx = self.derivatives[c]

                if dx > 0:
                    seq = self.contig.seq[pos:pos+30].replace("*", "")
                    side = "l"
                else:
                    seq = self.contig.seq[pos+1-30:pos+1].replace("*", "")
                    side = "r"
                print(f">{self.contig.name}_{c}_{pos}_{side}\n{seq}", file=fasta)


class SwitchpointFinder:
    def __init__(
        self, input_fn, outdir="./", window_size=7, min_depth=10, min_read_prop=0.01
    ):

        self.acefile = AceFile(input_fn)
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.outdir = outdir
        self.results = {}

    def fit(self):
        contig_dict = {}

        with open(f"{self.outdir}/boundaries_from_contigs.fas", "w") as self.fasta:
            all_reads = {}
            for _ in range(self.acefile.ncontigs):
                self.current_ctg = next(self.acefile)
                ctg = self.current_ctg
                for read in self.current_ctg.reads:
                    all_reads[read.name] = read
                if self.current_ctg.nreads / self.acefile.nreads > self.min_read_prop:
                    print(self.current_ctg.name)
                    cands, derivs = self.find_candidates()
                    contig_dict[self.current_ctg.name] = Result(ctg, cands, derivs)
                    self.generate_output_for_boundaries(cands, derivs)

                    for i in numpy.flatnonzero(cands):
                        d = derivs[i]
                        reads = get_reads_from_candidate(ctg, i)
                        for read in reads:
                            read.side = -numpy.sign(d).astype(int)
                            all_reads[read.name] = read
                        self.results[ctg.name] = SingleResult(ctg, i, d, ctg.depth[i], reads)

        for k in contig_dict:
            contig_dict[k].n = contig_dict[k].candidates.nonzero()[0].size

        return contig_dict, all_reads

    def generate_output_for_boundaries(self, cands, slopes):
        ctg = self.current_ctg
        contig_plot = ctg.generate_figure()
        for c, d in zip(numpy.flatnonzero(cands), slopes[cands]):
            boundary = Boundary(ctg, c, numpy.sign(d), abs(d))
            boundary.set_boundary_sequence()
            ctg.boundaries.append(boundary)
            print(boundary.boundary_seq_as_fasta(), file=self.fasta)
            logo_out = f"{self.outdir}/{ctg.name}_{boundary.side_as_l_or_r()}_logo.png"
            boundary.set_logo(save=logo_out)
            ctg.add_candidate_switchpoint_to_fig(contig_plot, c)
            contig_plot.savefig(f"{self.outdir}/{ctg.name}.png")
        pyplot.close(contig_plot)

    def find_candidates(self):
        contig = self.current_ctg
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
                wdiff = contig.unmasked[c-side:c+side+1] -\
                     contig.masked[c-side:c+side+1]
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
        candidates = sc * dmask * derivatives_mask

        return candidates.astype(bool), derivatives
