#! /usr/bin/env python3

import sys

import numpy
from numpy.lib.stride_tricks import as_strided






from ace import AceFile, Contig

acefile = AceFile(sys.argv[1])

for x in range(25):
    y = next(acefile)
