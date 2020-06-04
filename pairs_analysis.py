#! /usr/bin/env python3

import sys
import pickle

from difflib import get_close_matches

from kit.contig import Read, Pair





    def uncomplement_reads(self):
        if self.f.comp == 'C':
            self.f._comp_seq()
        if self.r.comp == 'C':
            self.r._comp_seq()


reads_pickle = sys.argv[1]
min_k = 8


with open(reads_pickle, 'rb') as pkl:
    mate_pairs = pickle.load(pkl)

pairs_list = []
for _, r1, r2 in mate_pairs:
    pair = Pair(f=r1, r=r2) if r1.name.endswith('f') else Pair(f=r2, r=r1)
    pair.uncomplement_reads()
    pairs_list.append(pair)

for pair in pairs_list:
    pair.set_reference()
    pair.get_kmers_in_reads()
    if pair.f and pair.r:
        print(pair)
        for k in pair.objective.kmers:
            print(k)
            pr_keys = list(pair.reference.kmers.keys())
            matches = get_close_matches(k, pr_keys, cutoff=6.5/8)
            print(k, pair.objective.kmers[k], " ".join([f"{x} {pair.reference.kmers[x]}" for x in matches]))
