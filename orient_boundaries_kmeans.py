#! /usr/bin/env python3

import numpy

from sklearn.cluster import KMeans

from ssw.lib import CSmithWaterman
from kit.utils import rc

aligner = CSmithWaterman()
aligner.set_alignment_params()
clusterer = KMeans(n_clusters=2)

fn = "/home/jamie/D/Local/191204_MITE_training_CLs/CCALI_CL0036/almitey"

with open(f"{fn}/boundaries_from_contigs.fas") as fasta:
    lines = [line.strip() for line in fasta]

# METHOD NOTES
# choose top boundary sequence
# align all against this sequence (perhaps also the revcomp)
# use resulting 2D alignment scores for kmeans clustering
