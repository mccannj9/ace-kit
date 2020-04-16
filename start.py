#! /usr/bin/env python3

import sys
import re

regex = re.compile(
    r"^CO CL(?P<cluster>\d+)Contig(?P<number>\d+)"
    r" (?P<length>\d+) (?P<nreads>\d+) \d+ [UC]$"
)

class AceFile(object):
    def __init__(self, filename):
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


class Contig(object):
    def __init__(self, lines):
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

        self.af = [
            line for line in self.lines if line.startswith('AF')
        ]
        self.rd = [
            line for line in self.lines if line.startswith('RD')
        ]


    @property
    def name(self):
        return f"CL{self.cluster}Contig{self.number}"

    def __repr__(self):
        return f"{self.name} with length {self.length} and {self.nreads} reads"

    def __len__(self):
        return self.length

    @property
    def padded_length(self):
        return self.length
    
    def unpadded_length(self):
        return len(self.seq.replace("*", ""))

    def padding(self):
        return [
            idx for idx, ltr in enumerate(self.seq) if ltr == "*"
        ]

ace = AceFile(sys.argv[1])

for x in range(10):
    y = next(ace)