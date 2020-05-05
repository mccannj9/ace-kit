#! /usr/bin/env python3

import os
import sys

from matplotlib import pyplot

from kit.finder import SwitchpointFinder

try:
    fn = sys.argv[2]
except IndexError:
    fn = "./"

finder = SwitchpointFinder(sys.argv[1], fn)
acefile = finder.acefile

results_dict = finder.fit()