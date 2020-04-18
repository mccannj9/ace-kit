#! /usr/bin/env python3

import sys
import re

from typing import List
from collections import namedtuple
from operator import attrgetter

import numpy
import matplotlib.pyplot as plt
plt.style.use('bmh')

StringVector = List[str]

Read = namedtuple('Read', ('name', 'length', 'start', 'f', 't'))

regex = re.compile(
    r"^CO CL(?P<cluster>\d+)Contig(?P<number>\d+)"
    r" (?P<length>\d+) (?P<nreads>\d+) \d+ [UC]$"
)


class AceFile(object):
    def __init__(self, filename:str):
        self.filename = filename
        self.file = open(filename)
        as_line = self.file.readline().strip().split()
        self.ncontigs, self.nreads = int(as_line[1]), int(as_line[2])
        
        # read until first CO line
        self.file.readline()
        self.buffer = [self.file.readline().strip()]

    def __next__(self):
        line = self.file.readline().strip()
        while not(line.startswith("CO")):
            self.buffer.append(line)
            line = self.file.readline().strip()
        output = self.buffer
        self.buffer = [line]
        return Contig(output)
    
    def readline(self):
        return self.file.readline()

    def close(self):
        self.file.close()


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

        self.reads = sorted(self.reads, key=attrgetter('start'))

        self.min = self.reads[0].start
        self.max = self.reads[-1].start + self.reads[-1].length - 1
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
    
    def plot_profile(self):
        fig, ax = plt.subplots(2, sharex=True, sharey=True)
        fig.suptitle(f"{self.name} masking pattern")
        ax[0].step(
            numpy.arange(self.min, self.max+1), self.unmasked, c='firebrick'
        )
        ax[0].set_ylabel(f'Unmasked Sites')
        ax[1].step(
            numpy.arange(self.min, self.max+1), self.masked, c='steelblue'
        )
        ax[1].set_xlabel(f'Position in {self.name} assembly')
        ax[1].set_ylabel(f'Masked Sites')
        return fig


    def __repr__(self):
        return f"id={self.name}: len={self.length}, nreads={self.nreads}"

    def __len__(self):
        return self.length

ace = AceFile(sys.argv[1])

for x in range(25):
    y = next(ace)
