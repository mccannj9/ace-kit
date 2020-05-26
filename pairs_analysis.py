#! /usr/bin/env python3

import sys
import pickle

from dataclasses import dataclass
from difflib import get_close_matches


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


reads_pickle = sys.argv[1]
min_k = 7


with open(reads_pickle, 'rb') as pkl:
    mate_pairs = pickle.load(pkl)

pairs_list = []
for _, r1, r2 in mate_pairs:
    pair = Pair(f=r1, r=r2) if r1.name.endswith('f') else Pair(f=r2, r=r1)
    pairs_list.append(pair)

for pair in pairs_list:
    pair.set_reference()
    pair.get_kmers_in_reads()
    if pair.f and pair.r:
        print(pair)
        for k in pair.objective.kmers:
            print(k)
            matches = get_close_matches(k, pair.reference.kmers.keys(), cutoff=0.75)
            print(k, " ".join([str(pair.reference.kmers[x]) for x in matches]))
