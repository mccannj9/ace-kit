#! /usr/bin/env python3

import sys

from ssw.lib import CSmithWaterman
from kit.utils import rc

aligner = CSmithWaterman()
# use default alignment params
aligner.set_alignment_params()

with open(f"{sys.argv[1]}/boundaries_from_contigs.fas") as fasta:
    lines = [line.strip() for line in fasta]

# METHOD NOTES
# choose top boundary sequence
# align all against this sequence (perhaps also the revcomp)
# use resulting 2D alignment scores for kmeans clustering

top_id = lines[0]
top = lines[1]

results = {}

for i, (_id, seq) in enumerate(zip(lines[::2], lines[1::2])):
    res_0 = aligner.align_sequence_pair(seq, top)
    res_1 = aligner.align_sequence_pair(
        "".join([rc[x] for x in seq[::-1]]), top
    )

    if res_0['nScore'] > res_1['nScore']:
        results[_id] = 1
    elif res_0['nScore'] < res_1['nScore']:
        results[_id] = 0
    else:
        # ambiguous result, should not happen
        results[_id] = None
