"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""

import cPickle as pickle
import numpy as np
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
    def __init__(self, nsamples):
        self.nsamples = nsamples

    def normal_deviates(self, mean, sdev):
        np.random.normal(mean, sdev, self.nsamples)

    def distance_samples(self):
        pass


class DpdfCollection(object):
    """
    Collection of all DPDFs and their properties. The class also contains
    methods for calculating Monte Carlo sampled properties.
    """
    ver = 'v210'

    def __init__(self, posteriors):
        self.bgps = read_cat('bgps_v210_evo').set_index(self.ver + 'cnum')
        self.dpdf_props = read_cat('bgps_v210_dpdf_props').set_index(self.ver +
                                                                     'cnum')
        self.posteriors = posteriors

    def construct_dpdfs(self):
        pass


