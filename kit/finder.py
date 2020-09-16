
import os
import itertools
from typing import List, Tuple, NewType
from statistics import mean

import numpy

from matplotlib import pyplot

from kit.ace import AceFile, NewAceFile
from kit.contig import Boundary
from kit.new_contig import NewContig, NewBoundary, sides
from kit.utils import rc

from ssw.lib import CSmithWaterman

BoundaryVec = List[Boundary]
Boundaries = NewType('Boundaries', List[NewBoundary])
Contigs = NewType('Contigs', List[NewContig])


class NewSwitchpointFinder(object):
    def __init__(
        self,
        input_fn: str,
        outdir: str = "./",
        window_size: int = 7,
        min_depth: int = 10,
        min_read_prop: float = 0.01,
        min_rate_change: float = 0.10,
        overwrite: bool = False
    ) -> None:

        self.acefile = NewAceFile(input_fn)
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.min_rate_change = 0.10
        self.outdir = outdir
        self.results = {}
        self.mean_contig_length = None
        self.min_contig_length = None

    def fit(self) -> Tuple[Contigs, Boundaries]:
        all_contigs = Contigs([])
        all_boundaries = Boundaries([])

        for _ in range(self.acefile.ncontigs):
            contig = next(self.acefile)
            if contig.nreads / self.acefile.nreads > self.min_read_prop:
                print(contig.name)
                self.find_candidates(contig)
                all_contigs += [contig]
                all_boundaries += contig.boundaries

        all_boundaries.sort(key=lambda x: x.rate, reverse=True)
        return all_contigs, all_boundaries

    def find_candidates(self, contig: NewContig) -> None:
        # here we find out when the sign of the difference between
        # unmasked and masked counts changes to find candidate boundaries
        sign_diff = numpy.sign(contig.masking_diff)
        sign_change = ((numpy.roll(sign_diff, 1) - sign_diff) != 0).astype(int)

        side = self.window_size // 2
        derivatives = numpy.zeros(shape=sign_diff.shape, dtype=float)
        for c in numpy.flatnonzero(sign_change):
            # this if statement disallows boundaries on the edges, because they
            # don't contain full windows (windows are center-anchored -> 00100)
            if (c-side < 0) or (c + side + 1 > contig.depth.size):
                sign_change[c] = 0
            else:
                win_depth = contig.depth[c-side:c+side+1].mean()
                derivatives[c] = (
                    contig.masking_diff[c+side] - contig.masking_diff[c-side]
                ) / self.window_size

                if abs(derivatives[c]) / win_depth < self.min_rate_change:
                    sign_change[c] = 0

        # get the idx of min and max of derivative
        deriv_mask = numpy.zeros(shape=derivatives.shape)
        deriv_mask[derivatives.argmin()] = 1
        deriv_mask[derivatives.argmax()] = 1
        cands = sign_change * contig.depth_mask(self.min_depth) * deriv_mask
        contig.generate_figure()

        for c, d in zip(numpy.flatnonzero(cands), derivatives[cands.astype(bool)]):
            boundary = NewBoundary(
                contig, c, int(numpy.sign(d)), abs(d)
            )
            boundary.set_logo(
                save=f"{self.outdir}/{contig.name}_{boundary.side_as_l_or_r()}.png"
            )
            boundary.set_boundary_sequence()
            boundary.add_boundary_to_contig_profile_plot(contig)
            contig.boundaries.append(boundary)

        contig.fig.savefig(
            f"{self.outdir}/{contig.name}.png", figsize=(6, 4)
        )
        pyplot.close(contig.fig)

    def orient_boundaries(self, boundaries: Boundaries) -> List[str]:
        fasta_output = ""
        if len(boundaries) < 2:
            return fasta_output

        for b in boundaries:
            if b.side == 1:
                seq = b.seq
            elif b.side == -1:
                seq = "".join([rc[x] for x in b.seq[::-1]])

            fasta_output += f">{b.contig.name}_{sides[b.side]}\n{seq}\n"

        return fasta_output


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
        # derivatives_mask[min_der_idx] = 1
        # derivatives_mask[max_der_idx] = 1
        derivatives_mask[derivatives.argmin()] = 1
        derivatives_mask[derivatives.argmax()] = 1

        candidates = sc * dmask * derivatives_mask

        return candidates.astype(bool), derivatives

    def orient_boundaries(self, boundaries: BoundaryVec, **params):
        if len(boundaries) < 1:
            return []
        aligner = CSmithWaterman()
        aligner.set_alignment_params(**params)

        # assume boundaries are sorted
        top = boundaries[0]
        oriented = []

        for b in boundaries[1:]:
            res_0 = aligner.align_sequence_pair(b.seq, top.seq)
            res_1 = aligner.align_sequence_pair(
                "".join([rc[x] for x in b.seq[::-1]]), top.seq
            )
            print(res_0['nScore'], res_1['nScore'], b.contig_name)

            if res_0['nScore'] > res_1['nScore']:
                oriented.append(b.seq)
                b.orient = 0
            elif res_0['nScore'] < res_1['nScore']:
                oriented.append(
                    "".join([rc[x] for x in b.seq[::-1]])
                )
                b.orient = 1
            else:
                # ambiguous result, should not happen
                b.orient = -999

        return oriented

    def estimate_TIR_length(self, boundaries: BoundaryVec, **align_params):
        orientations = set([
            b.orient for b in boundaries
        ])

        if len(boundaries) < 2 and len(orientations) != 2:
            print("Cannot estimate TIR length")
            return None

        ext = self.min_contig_length // 2
        ext = 80

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
            l = abs(results['nRefEnd'] - results['nRefBeg'])
            print(x.contig, x.side, x.orient, y.contig, y.side, y.orient, l)
            alignment_lengths.append(
                abs(results['nRefEnd'] - results['nRefBeg'])
            )
            print(results['sAlignment'])
        return mean(alignment_lengths), alignment_lengths
