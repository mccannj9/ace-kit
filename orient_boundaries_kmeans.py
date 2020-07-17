#! /usr/bin/env python3

import sys

import numpy

from sklearn.cluster import KMeans

from ssw.lib import CSmithWaterman
from kit.utils import rc

aligner = CSmithWaterman()
aligner.set_alignment_params()
clusterer = KMeans(n_clusters=2)

with open(f"{sys.argv[1]}/boundaries_from_contigs.fas") as fasta:
    lines = [line.strip() for line in fasta]

# METHOD NOTES
# choose top boundary sequence
# align all against this sequence (perhaps also the revcomp)
# use resulting 2D alignment scores for kmeans clustering

top = lines[1]
data = numpy.zeros(shape=(len(lines)//2, 2), dtype=numpy.float)


for i, seq in enumerate(lines[1::2]):
    data[i, 0] = aligner.align_sequence_pair(seq, top)['nScore']
    data[i, 1] = aligner.align_sequence_pair(
        "".join([rc[x] for x in seq[::-1]]), top
    )['nScore']

preds = clusterer.fit_predict(data)