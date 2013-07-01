"""
==========================
Kinetic Temperatures Plots
==========================

Plots for kinetic temperatures (T_K) and other line properties derived from
NH3 observations.

"""
# TODO add plots for line widths and intensities

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from .. import catalog
from .. import util
import ipdb as pdb

def ktemp_hist(bgps=[]):
