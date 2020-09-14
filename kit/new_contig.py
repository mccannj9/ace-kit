
import re

from dataclasses import dataclass
from typing import List
from operator import attrgetter

import numpy
import logomaker

from matplotlib import pyplot
from logomaker import Logo

from kit.utils import compute_end_pos, window, create_seqlogo_dataframe, colors
from kit.utils import sign
from kit.html import minor_row_template

regex = re.compile(
    r"^CO CL(?P<cluster>\d+)Contig(?P<number>\d+)"
    r" (?P<length>\d+) (?P<nreads>\d+) \d+ [UC]$"
)


pyplot.style.use('bmh')

StringVector = List[str]
sides = ["X", "left", "right"]

@dataclass
class NewRead:
    name: str
    length: int
    start: str
    f: int
    t: int
    seq: str
    ctg: str
    comp: str

    @property
    def end(self):
        return self.start + self.length

    @property
    def extent(self):
        return range(self.start, self.end)

    def overlapping(self, pos):
        return pos in self.extent

    def __len__(self):
        return len(self.seq)


ReadsVec = List[NewRead]


def get_sequences_for_contig(ctg):
    seqs = []
    for i, line in enumerate(ctg.lines):
        if line.startswith('RD'):
            seq = []
            count = 1
            cl = ctg.lines[i+count]
            while cl:
                seq.append(cl)
                count += 1
                cl = ctg.lines[i+count]
            seqs.append("".join(seq))
    return seqs


class NewContig(object):
    def __init__(self, lines: StringVector) -> None:
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
            NewRead(
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
            r.start += self.shift
            unmasked = numpy.zeros(r.length).astype(bool)
            unmasked[r.f-1:r.t-1] = True
            masked = ~unmasked
            self.unmasked[r.start: r.end] += unmasked
            self.masked[r.start:r.end] += masked

        self.masking_diff = self.unmasked - self.masked

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

    def depth_mask(self, min_depth):
        return self.depth >= min_depth

    def plot_profile(self, ax, window_size=1, shift=1):
        if window_size > 1:
            unmasked_data = window(self.unmasked, window=window_size, shift=shift).mean(axis=1)
            masked_data = window(self.masked, window=window_size, shift=shift).mean(axis=1)
        else:
            unmasked_data = self.unmasked
            masked_data = self.masked

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
            xaxis, self.depth, c='seagreen', label='Depth'
        )
        ax.legend()
        return artist1, artist2, artist3

    def generate_figure(self, fs=(6, 4), filename=None, **kwargs):
        self.fig, ax = pyplot.subplots(1, figsize=fs)
        self.plot_profile(ax)
        self.fig.suptitle(f"{self.name} RD = {round(self.average_rd, 1)}")
        if filename:
            self.fig.savefig(filename, **kwargs)
            pyplot.close(self.fig)

    def add_candidate_switchpoint_to_fig(self, fig, candidate):
        ymax = self.depth.max()
        for ax in fig.axes:
            ax.vlines(self.min + candidate, 0, ymax, linestyles='dotted')

    def __repr__(self):
        return f"id={self.name}: len={self.length}, nreads={self.nreads}"

    def __len__(self):
        return self.length

    def reads_on_position(self, pos: int) -> ReadsVec:
        overlapping_reads = []
        for r in self.reads:
            if r.overlapping(pos):
                overlapping_reads.append(r)
        return overlapping_reads


class NewBoundary:
    def __init__(
        self,
        contig: NewContig,
        pos: int,
        side: int,
        rate: float,
    ) -> None:
        self.contig = contig
        self.contig_name = contig.name
        self.pos = pos
        self.depth = self.contig.depth[self.pos]
        self.side = side
        self.rate = rate

        # unset attributes set by instance methods
        self.reads = None
        self.logo = None
        self.orient = None
        self.seq = None

        self.add_boundary_to_contig_profile_plot(contig)

    @property
    def name(self):
        return f"{self.contig.name}_{self.side_as_l_or_r()}"

    @property
    def slope(self):
        return self.rate * self.side

    def set_boundary_sequence(self, after: int=30) -> None:
        pos = self.pos - self.contig.shift
        if self.side == 1:
            self.seq = self.contig.seq[pos-1:pos+after-1].replace("*", "")
            self.seq
        else:
            self.seq = self.contig.seq[pos-after:pos].replace("*", "")

    def set_overlapping_reads(self) -> None:
        self.reads = self.contig.reads_on_position(self.pos)

    def get_boundary_seqs_from_reads(self, before: int=5, after: int=30) -> StringVector:
        seqs = []

        if not(self.reads):
            self.set_overlapping_reads()

        for read in self.reads:
            pos_in_read = self.pos - read.start
            if self.slope > 0:
                left = pos_in_read - before
                right = pos_in_read + after
            else:
                left = pos_in_read - after + 1
                right = pos_in_read + before + 1

            # have to check if we need to add a padding char
            # then extracted seqs should all have same length
            if left < 0:
                pad_left = "0" * abs(left)
                left = 0

            else:
                pad_left = ""

            if right > (read.length - 1):
                pad_right = "0" * abs(right - (read.length))
                right = read.length

            else:
                pad_right = ""

            seqs.append(
                f"{pad_left}{read.seq[left:right]}{pad_right}".replace("*", "-")
            )

        return seqs

    def add_boundary_to_contig_profile_plot(self, contig: NewContig) -> None:
        ymax = contig.depth.max()
        for ax in contig.fig.axes:
            ax.vlines(self.pos - self.contig.shift, 0, ymax, linestyles='dotted')

    def set_logo(
        self,
        before: int = 5,
        after: int = 30,
        figsize=(10, 2.5),
        save=None
    ) -> None:
        seqs = self.get_boundary_seqs_from_reads()
        counts, freqs, df = create_seqlogo_dataframe(seqs)

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

    def table_row_template(self):
        d = self.__dict__.copy()
        d['side'] = sides[int(self.side)]
        d['pos'] = self.pos - self.contig.shift
        d['rate'] = round(d['rate'], 1)
        return minor_row_template.safe_substitute(d)
