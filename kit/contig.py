
import re

from typing import List
from collections import namedtuple
from operator import attrgetter

import numpy

from matplotlib import pyplot

from kit.utils import compute_end_pos, window

pyplot.style.use('bmh')

StringVector = List[str]

Read = namedtuple('Read', ('name', 'length', 'start', 'f', 't'))

regex = re.compile(
    r"^CO CL(?P<cluster>\d+)Contig(?P<number>\d+)"
    r" (?P<length>\d+) (?P<nreads>\d+) \d+ [UC]$"
)


class Contig(object):
    def __init__(self, lines:StringVector):
        # just initializing values so editor doesn't complain
        self.length = 0
        self.nreads = 0
        self.cluster = 0
        self.number = 0

        intro = lines[0]
        self.lines = lines
        
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

        lengths = [
            int(line.split()[2]) for line in self.lines if line.startswith('RD')
        ]

        froms = [
            int(line.split()[1]) for line in self.lines if line.startswith('QA')
        ]

        tos = [
            int(line.split()[2]) for line in self.lines if line.startswith('QA')
        ]

        zipper = zip(names, lengths, starts, froms, tos)

        self.reads = [
            Read(
                name=_id, length=l, start=s, f=f, t=t
            ) for _id, l, s, f, t in zipper
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
        depth =  unmasked_data + masked_data

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
        # self.plot_profile(ax[1], window_size=3)
        # self.plot_profile(ax[2], window_size=5)
        fig.suptitle(f"{self.name} RD = {round(self.average_rd, 1)}")
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        if filename:
            fig.savefig(filename, **kwargs)
        return fig
    
    def add_candidate_switchpoints_to_fig(self, fig, candidates):
        b, e = candidates
        # gets max depth for drawing the line
        ymax = (self.unmasked + self.masked).max()
        for i, ax in enumerate(fig.axes):
            ax.vlines(self.min + b, 0, ymax, linestyles='dotted')
            ax.vlines(self.max + e - (2*i + 1), 0, ymax, linestyles='dotted')
        

    def __repr__(self):
        return f"id={self.name}: len={self.length}, nreads={self.nreads}"

    def __len__(self):
        return self.length