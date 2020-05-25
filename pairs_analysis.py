#! /usr/bin/env python3

import sys
from dataclasses import dataclass
from difflib import get_close_matches

@dataclass
class Read:
    name: str
    seq: str
    boundary: int = 0

    def __bool__(self):
        return self.boundary == 1

    def __mul__(self, other):
        return self.boundary * other.boundary

    def __len__(self):
        return len(self.seq)

@dataclass
class Pair:
    f: Read
    r: Read

    def __bool__(self):
        return self.f * self.r

    def set_reference(self):
        self.reference = self.f if self.f else self.r
        self.objective = self.f if self.reference != self.f else self.r

    def get_kmers_in_reads(self, k=7):
        self.f.kmers = {self.f.seq[x:x+k]:x for x in range(len(self.f)-k+1)}
        self.r.kmers = {self.r.seq[x:x+k]:x for x in range(len(self.r)-k+1)}


def create_pair_from_list(seqlist):
    names = seqlist[0], seqlist[2]
    id1, id2 = names[0][1:].split("_")[0], names[1][1:].split("_")[0]
    bd1, bd2 = names[0].split("_")[-1], names[1].split("_")[-1]
    r1 = Read(name=id1, seq=seqlist[1], boundary=int(bd1))
    r2 = Read(name=id2, seq=seqlist[3], boundary=int(bd2))
    return Pair(f=r1, r=r2) if r1.name.endswith('f') else Pair(f=r2, r=r1)


reads_fn = sys.argv[1]
min_k = 7

pairs_list = []
with open(reads_fn) as fas:
    seqs = []
    for i, line in enumerate(fas, start=1):
        seqs.append(line.rstrip())
        if i % 4 == 0:
            pairs_list.append(create_pair_from_list(seqs))
            seqs = []

for pair in pairs_list:
    pair.set_reference()
    pair.get_kmers_in_reads()

for k in pair.objective.kmers:
    print(k)
    matches = get_close_matches(k, pair.reference.kmers.keys(), cutoff=0.75)
    print(k, " ".join([str(pair.reference.kmers[x]) for x in matches]))
