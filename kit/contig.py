
import os
import re

from dataclasses import dataclass
from typing import List
from operator import attrgetter

import numpy
import logomaker

from matplotlib import pyplot
from logomaker import Logo

from kit.utils import compute_end_pos, window, rc, create_seqlogo_dataframe, colors
from kit.html import minor_row_template


pyplot.style.use('bmh')

StringVector = List[str]


@dataclass
class Read:
    name: str
    length: int
    start: str
    f: int
    t: int
    seq: str
    ctg: str
    comp: str
    boundary: int = 0
    side: int = 0

    def _fasta(self):
        seq = self.seq.replace("*", "")
        if self.comp == 'C':
            self._comp_seq()
        return f">{self.name}_{self.f}_{self.t}_{self.boundary}\n{seq}"

    def _comp_seq(self):
        self.seq = "".join([rc[c] for c in self.seq[::-1]])
        self.comp = 'U' if self.comp == 'C' else 'C'
        if self.boundary:
            self.boundary = len(self.seq) - self.boundary - 1

    def _boundary_seq(self, length=30):
        if self.side == 0:
            return self.seq
        return self.seq[self.boundary:self.boundary+length]

    def __bool__(self):
        return self.boundary != 0

    def __mul__(self, other):
        return self.boundary * other.boundary

    def __len__(self):
        return len(self.seq)


ReadsVec = List[Read]


@dataclass
class Pair:
    f: Read
    r: Read

    def __bool__(self):
        return self.f * self.r != 0

    def set_reference(self):
        self.reference = self.f if self.f else self.r
        self.objective = self.f if self.reference != self.f else self.r

    def get_kmers_in_reads(self, k=7, length=30):
        # bseq_f = self.f.seq[self.f.boundary]
        bseq_f = self.f._boundary_seq(length=length)
        bseq_r = self.r._boundary_seq(length=length)
        self.f.kmers = {bseq_f[x:x+k]:x for x in range(len(bseq_f)-k+1)}
        self.r.kmers = {bseq_r[x:x+k]:x for x in range(len(bseq_r)-k+1)}

    def uncomplement_reads(self):
        if self.f.comp == 'C':
            self.f._comp_seq()
        if self.r.comp == 'C':
            self.r._comp_seq()


regex = re.compile(
    r"^CO CL(?P<cluster>\d+)Contig(?P<number>\d+)"
    r" (?P<length>\d+) (?P<nreads>\d+) \d+ [UC]$"
)


def get_sequences_for_contig(ctg):
    seqs = []
    for i, line in enumerate(ctg.lines):
        if line.startswith('RD'):
            seq = []
            count = 1
            l = ctg.lines[i+count]
            # while l:= ctg.lines[i+count]:
            while l:
                seq.append(l)
                count += 1
                l = ctg.lines[i+count]
            seqs.append("".join(seq))
    return seqs


class Contig(object):
    def __init__(self, lines: StringVector):
        # just initializing values so editor doesn't complain
        self.length = 0
        self.nreads = 0
        self.cluster = 0
        self.number = 0

        intro = lines[0]
        self.lines = lines
        self.boundaries = []

        matches = regex.match(intro).groupdict()
        for m in matches:
            setattr(self, m, int(matches[m]))

        first_empty = self.lines.index('')
        self.seq = "".join(self.lines[1:first_empty])
        self.lines = self.lines[first_empty+1:]

        first_empty = self.lines.index('')
        self.bq = " ".join(self.lines[1:first_empty]).split()
        # pad the bq with -1 where seq = *
        for x in self.padding():
            self.bq.insert(x, -1)

        self.lines = self.lines[first_empty+1:]

        names = [
            line.split()[1] for line in self.lines if line.startswith('AF')
        ]

        starts = [
            int(line.split()[3]) for line in self.lines if line.startswith('AF')
        ]

        comps = [
            line.split()[2] for line in self.lines if line.startswith('AF')
        ]

        lengths = [
            int(line.split()[2]) for line in self.lines if line.startswith('RD')
        ]

        froms = [
            int(line.split()[1]) for line in self.lines if line.startswith('QA')
        ]

        tos = [
            int(line.split()[2]) for line in self.lines if line.startswith('QA')
        ]

        seqs = get_sequences_for_contig(self)

        zipper = zip(names, lengths, starts, froms, tos, seqs, comps)

        self.reads = [
            Read(
                name=_id, length=l, start=s, f=f, t=t, seq=_seq, ctg=self.name, comp=c
            ) for _id, l, s, f, t, _seq, c in zipper
        ]

        # have to sort by starting and ending position
        # sometimes (with gaps) sorting gives different results
        reads_for_min = sorted(self.reads, key=attrgetter('start'))
        reads_for_max = sorted(self.reads, key=compute_end_pos)

        self.min = reads_for_min[0].start
        self.max = reads_for_max[-1].start + reads_for_max[-1].length - 1
        self.shift = abs(self.min) if self.min < 1 else -1
        self.assembly_len = self.shift + self.max + 1 if self.min < 1 else self.max

        # sum over read positions
        # use from-to, start and shift to calculate an array of zeros
        # and ones for positions in the assembly
        # add this to unmasked, invert and add inversion to masked
        self.unmasked = numpy.zeros(shape=self.assembly_len, dtype=numpy.int64)
        self.masked = numpy.zeros(shape=self.assembly_len, dtype=numpy.int64)
        for r in self.reads:
            unmasked = numpy.zeros(r.length).astype(bool)
            unmasked[r.f-1:r.t-1] = True
            masked = ~unmasked
            begin = self.shift + r.start
            end = self.shift + r.start + r.length
            self.unmasked[begin: end] += unmasked
            self.masked[begin:end] += masked

        self.depth = self.unmasked + self.masked
        self.average_rd = self.depth.mean()

    @property
    def nboundaries(self):
        return len(self.boundaries)

    @property
    def boundary_rate_sum(self):
        return sum([x.rate for x in self.boundaries])

    @property
    def name(self):
        return f"CL{self.cluster}Contig{self.number}"

    @property
    def padded_length(self):
        return self.length

    @property
    def unpadded_length(self):
        return len(self.seq.replace("*", ""))

    def padding(self):
        return [
            idx for idx, ltr in enumerate(self.seq) if ltr == "*"
        ]

    def plot_profile(self, ax, window_size=1, shift=1):
        if window_size > 1:
            unmasked_data = window(self.unmasked, window=window_size, shift=shift).mean(axis=1)
            masked_data = window(self.masked, window=window_size, shift=shift).mean(axis=1)
        else:
            unmasked_data = self.unmasked
            masked_data = self.masked
        depth = unmasked_data + masked_data

        ax.set_ylabel('No. of Sites')

        ax.set_xlabel(f'Position, Window Size = {window_size}')
        xaxis = numpy.arange(self.min, self.max-window_size+2)
        artist1 = ax.step(
            xaxis, unmasked_data, c='firebrick', label='Unmasked'
        )
        artist2 = ax.step(
            xaxis, masked_data, c='steelblue', label='Masked'
        )
        artist3 = ax.step(
            xaxis, depth, c='seagreen', label='Depth'
        )
        ax.legend()
        return artist1, artist2, artist3

    def generate_figure(self, fs=(6, 4), filename=None, **kwargs):
        fig, ax = pyplot.subplots(1, figsize=fs)
        self.plot_profile(ax)
        fig.suptitle(f"{self.name} RD = {round(self.average_rd, 1)}")
        if filename:
            fig.savefig(filename, **kwargs)
            pyplot.close(fig)
        return fig

    def add_candidate_switchpoint_to_fig(self, fig, candidate):

        ymax = self.depth.max()
        for ax in fig.axes:
            ax.vlines(self.min + candidate, 0, ymax, linestyles='dotted')

    def __repr__(self):
        return f"id={self.name}: len={self.length}, nreads={self.nreads}"

    def __len__(self):
        return self.length

    def get_reads_with_position(self, pos, slope, l=30):
        seqs = []
        for r in self.reads:
            b = r.start
            e = r.start + r.length
            if ((pos - self.shift) in range(b, e)):
                p = pos - self.shift - b
                seq = r.seq.replace("*", "-")
                if slope > 0:
                    seqs.append(seq[p:p+l])
                else:
                    seqs.append(seq[p-l+1:p+1])
        # this removes those reads which do not have length == l
        # otherwise dataframe is messed up
        return [x for x in seqs if len(x) == l]


@dataclass
class Boundary:
    contig: Contig
    contig_name: str
    pos: int
    depth: int
    side: int
    rate: float
    logo: Logo = None
    orient: int = 0
    seq: str = ""

    @property
    def name(self):
        return f"{self.contig.name}_{self.side_as_l_or_r()}"

    def set_boundary_sequence(self, l=30):
        pos = self.pos - self.contig.shift
        if self.side == 1:
            self.seq = self.contig.seq[pos:pos+l].replace("*", "")
        else:
            self.seq = self.contig.seq[pos+1-l:pos+1].replace("*", "")

    def set_logo(self, l=30, figsize=(10, 2.5), save=None):
        seqs = self.contig.get_reads_with_position(self.pos, self.side * self.rate, l=30)
        df = create_seqlogo_dataframe(seqs)

        fig, ax = pyplot.subplots(1, figsize=figsize)
        self.logo = logomaker.Logo(df, color_scheme=colors, ax=ax)
        self.logo.style_xticks(anchor=0, spacing=5, rotation=45)
        self.logo.ax.set_xlim([-1, len(df)])
        self.logo.ax.set_ylabel('information (bits)')

        if save:
            fig.savefig(save)
            self.logo_path = save
            pyplot.close(fig)

    def side_as_l_or_r(self):
        return "l" if self.side == 1 else "r"

    def boundary_seq_as_fasta(self):
        return f">{self.contig.name}_{self.side_as_l_or_r()}\n{self.seq}"

    def get_reads_from_boundary(self, buffer=0):
        reads = []
        for r in self.contig.reads:
            b = r.start
            e = r.start + r.length
            # read_range = range(b, e)
            # true_pos = self.pos - self.contig.shift
            # contig_range = range(true_pos, true_pos + buffer)
            # overlap = range(
            #     max(read_range[0], contig_range[0]), max(read_range[-1], contig_range[-1])
            # )
            # if overlap.start < overlap.stop:
            if ((self.pos - self.contig.shift) in range(b, e)):
                reads.append(r)
        return reads

    def get_mate_pairs(self, reads: ReadsVec):
        all_reads = {x.name: x for x in self.contig.reads}
        pairs_dict = {}

        for r in reads:
            r_end = r.name[-1]
            mate_end = "f" if r_end == "f" else "r"
            mate_id = f"{r.name[:-1]}{mate_end}"
            if mate_id in all_reads:
                if r.name[:-1] not in pairs_dict:
                    pairs_dict[r.name[:-1]] = [r, all_reads[mate_id]]
                    pairs_dict[r.name[:-1]].sort(key=lambda x: x.name[:-1])

    def table_row_template(self):
        d = self.__dict__
        d['side'] = self.side_as_l_or_r()
        d['side'] = "left" if d['side'] == "l" else "right"
        d['pos'] = self.pos - self.contig.shift
        d['rate'] = round(d['rate'], 1)
        return minor_row_template.safe_substitute(d)
