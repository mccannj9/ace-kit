
import numpy

from matplotlib import pyplot

from kit.ace import AceFile
from kit.contig import Boundary
from kit.utils import rc

from ssw.lib import CSmithWaterman


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
        contigs_list = []

        with open(f"{self.outdir}/boundaries_from_contigs.fas", "w") as self.fasta:
            all_reads = {}
            for _ in range(self.acefile.ncontigs):
                self.current_ctg = next(self.acefile)
                all_reads.update({x.name: x for x in self.current_ctg.reads})
                if self.current_ctg.nreads / self.acefile.nreads > self.min_read_prop:
                    print(self.current_ctg.name)
                    cands, derivs = self.find_candidates()
                    self.generate_output_for_boundaries(cands, derivs)
                    contigs_list.append(self.current_ctg)

        return contigs_list, all_reads

    def generate_output_for_boundaries(self, cands, slopes):
        ctg = self.current_ctg
        contig_plot = ctg.generate_figure()
        for c, d in zip(numpy.flatnonzero(cands), slopes[cands]):

            boundary = Boundary(ctg, ctg.name, c, ctg.depth[c], numpy.sign(d), abs(d))
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

    def orient_boundaries(self, boundaries, **params):
        aligner = CSmithWaterman()
        aligner.set_alignment_params(**params)
        boundaries.sort(key=lambda x: x.rate, reverse=True)

        top = boundaries[0].seq
        ids = [b._id for b in boundaries]
        seqs = [b.seq for b in boundaries]
        results = {}

        for _id, seq in zip(ids, seqs):
            res_0 = aligner.align_sequence_pair(seq, top)
            res_1 = aligner.align_sequence_pair(
                "".join([rc[x] for x in seq[::-1]]), top
            )

            if res_0['nScore'] > res_1['nScore']:
                results[_id] = 1
            elif res_0['nScore'] < res_1['nScore']:
                results[_id] = 0
            else:
                # ambiguous result, should not happen
                results[_id] = None
