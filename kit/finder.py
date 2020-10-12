
import os
import itertools
from typing import List, Tuple, NewType, TextIO
from statistics import mean

import numpy

from matplotlib import pyplot

from kit.ace import AceFile
from kit.contig import Boundary, Contig, sides
from kit.utils import rc, revcomp, unique_kmer_distance
from kit.utils import first_nonzero, last_nonzero

from ssw.lib import CSmithWaterman, default_alignment_parameters

Boundaries = NewType('Boundaries', List[Boundary])
Contigs = NewType('Contigs', List[Contig])


class SwitchpointFinder(object):
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

        self.acefile = AceFile(input_fn)
        self.window_size = window_size
        self.min_depth = min_depth
        self.min_read_prop = min_read_prop
        self.min_rate_change = 0.10
        self.outdir = outdir
        self.results = {}
        self.mean_contig_length = None
        self.min_contig_length = None

    def fit(self, log: TextIO) -> Tuple[Contigs, Boundaries]:
        all_contigs = Contigs([])
        all_boundaries = Boundaries([])

        for _ in range(self.acefile.ncontigs):
            contig = next(self.acefile)
            if contig.nreads / self.acefile.nreads > self.min_read_prop:
                print(contig.name)
                self.find_candidates(contig)
                all_contigs += [contig]
                all_boundaries += contig.boundaries

        tirs = self.detect_tirs(all_contigs, **default_alignment_parameters)

        for ctg, res in tirs:
            print(ctg.name, ctg.length, res.best_score, res.alignment_length, file=log)
            print(res.alignment_repr + "\n", file=log)

        all_boundaries.sort(key=lambda x: x.rate, reverse=True)
        return all_contigs, all_boundaries

    def find_candidates(self, contig: Contig) -> None:
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
        # enforces criterium that min must be negative, mask diff decreasing (right bound)
        deriv_mask[derivatives.argmin()] = 1 if derivatives[derivatives.argmin()] < 0 else 0
        # enforces criterium that max must be positive, mask diff increasing (left bound)
        deriv_mask[derivatives.argmax()] = 1 if derivatives[derivatives.argmax()] > 0 else 0

        new_derivatives = derivatives * deriv_mask
        # disallow min derivative to be to the left of max derivative
        # left bound must come before right bound, obviously
        if new_derivatives.argmin() < new_derivatives.argmax():
            x = abs(new_derivatives)
            idx = numpy.where(x != 0, x, x[x.argmax()]).argmin()
            deriv_mask[idx] = 0

        cands = sign_change * contig.depth_mask(self.min_depth) * deriv_mask
        contig.generate_figure()

        for c, d in zip(numpy.flatnonzero(cands), derivatives[cands.astype(bool)]):
            boundary = Boundary(
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

    def detect_tirs(self, contigs: Contigs, **align_params):
        results = []

        for contig in contigs:
            if len(contig.boundaries) == 2:
                print(f"Checking {contig.name} for TIR")
                aligner = CSmithWaterman(debug=False)
                aligner.set_alignment_params(**align_params)
                ssw_result = aligner.align_sequence_pair(
                    contig.boundaries[0].seq, revcomp(contig.boundaries[-1].seq)
                )
                print(unique_kmer_distance(
                    contig.boundaries[0].seq, revcomp(contig.boundaries[-1].seq), 5
                ))
                print(ssw_result.alignment_repr)

                if ssw_result.alignment_length >= 10:
                    results.append((contig, ssw_result))

        return results


    def estimate_TIR_length(self, boundaries: Boundaries, **align_params):
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
