
import itertools
from typing import List
from statistics import mean

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
        self.mean_contig_length = None
        self.min_contig_length = None

    def fit(self):
        contigs_list = []

        with open(f"{self.outdir}/boundaries_from_contigs.fas", "w") as self.fasta:
            all_reads = {}
            all_boundaries = []
            for _ in range(self.acefile.ncontigs):
                self.current_ctg = next(self.acefile)
                all_reads.update({x.name: x for x in self.current_ctg.reads})
                if self.current_ctg.nreads / self.acefile.nreads > self.min_read_prop:
                    print(self.current_ctg.name)
                    cands, derivs = self.find_candidates()
                    boundaries = self.generate_output_for_boundaries(cands, derivs)
                    contigs_list.append(self.current_ctg)
                    all_boundaries += boundaries
            all_boundaries.sort(key=lambda x: x.rate, reverse=True)
            oriented_seqs = self.orient_boundaries(all_boundaries)

        contig_lengths = [c.length for c in contigs_list if c.boundaries]

        if contig_lengths:
            self.mean_contig_length = mean(contig_lengths)
            self.min_contig_length = min(contig_lengths)

        return contigs_list, all_boundaries, oriented_seqs, all_reads

    def generate_output_for_boundaries(self, cands, slopes):
        ctg = self.current_ctg
        contig_plot = ctg.generate_figure()
        boundaries = []
        for c, d in zip(numpy.flatnonzero(cands), slopes[cands]):

            boundary = Boundary(ctg, ctg.name, c, ctg.depth[c], numpy.sign(d), abs(d))
            boundary.set_boundary_sequence()
            ctg.boundaries.append(boundary)
            print(boundary.boundary_seq_as_fasta(), file=self.fasta)
            logo_out = f"{self.outdir}/{ctg.name}_{boundary.side_as_l_or_r()}_logo.png"
            boundary.set_logo(save=logo_out)
            ctg.add_candidate_switchpoint_to_fig(contig_plot, c)
            contig_plot.savefig(f"{self.outdir}/{ctg.name}.png")
            boundaries.append(boundary)
        pyplot.close(contig_plot)

        return boundaries

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

    def orient_boundaries(self, boundaries:BoundaryVec, **params):
        if len(boundaries) < 1:
            return []
        aligner = CSmithWaterman()
        aligner.set_alignment_params(**params)

        # assume boundaries are sorted
        top = boundaries[0]
        oriented = []
        results = {}

        for b in boundaries[1:]:
            res_0 = aligner.align_sequence_pair(b.seq, top.seq)
            res_1 = aligner.align_sequence_pair(
                "".join([rc[x] for x in b.seq[::-1]]), top.seq
            )

            if res_0['nScore'] > res_1['nScore']:
                results[b._id] = 1
                oriented.append(b.seq)
                b.orient = 0
            elif res_0['nScore'] < res_1['nScore']:
                results[b._id] = 0
                oriented.append(
                    "".join([rc[x] for x in b.seq[::-1]])
                )
                b.orient = 1
            else:
                # ambiguous result, should not happen
                results[b._id] = None
                b.orient = -999

        return oriented

    def estimate_TIR_length(self, boundaries:BoundaryVec, **align_params):
        orientations = set([
            b.orient for b in boundaries
        ])

        if len(boundaries) < 2 and len(orientations) != 2:
            print("Cannot estimate TIR length")
            return None

        ext = self.min_contig_length // 2

        aligner = CSmithWaterman(debug=False)
        aligner.set_alignment_params(**align_params)

        b0 = [b for b in boundaries if not(b.orient)]
        b1 = [b for b in boundaries if b.orient]
        alignment_lengths = []

        for x, y in itertools.product(b0, b1):
            seq1 = x.extract_seq_from_contig(ext)
            # extract and revcomp
            seq2 = y.extract_seq_from_contig(ext)
            seq2 = "".join([rc[x] for x in seq2[::-1]])
            results = aligner.align_sequence_pair(seq1, seq2)
            print(x.contig, x.side, x.orient, y.contig, y.side, y.orient)
            alignment_lengths.append(
                abs(results['nRefEnd'] - results['nRefBeg'])
            )
            print(results['sAlignment'])
        return mean(alignment_lengths), alignment_lengths
