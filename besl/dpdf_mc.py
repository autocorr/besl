"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""

import cPickle as pickle
from .catalog import read_cat, read_dpdf
from .coord import sep


class Dpdf(object):
    """
    Container object for a BGPS source with a DPDF
    """
    def __init__(self):
        # pdf
        # kinetic temperature
        # initiate sampler
        pass


class Sampler(object):
    # x-dist values
    def __init__(self):
        pass


class DpdfCollection(object):
    """
    Collection of all DPDFs and their properties. The class also contains
    methods for calculating Monte Carlo sampled properties.
    """
    bgps = None
    dpdf_cat = None
    ppv_group = None
    dpdfs = {}

    def __init__(self):
        pass

    def construct_dpdfs(self):
        pass



