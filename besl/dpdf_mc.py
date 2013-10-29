"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""

import numpy as np
from .catalog import read_cat, read_dpdf
from scipy.stats import gaussian_kde


class ClumpProp(object):
    """
    Container object for a BGPS source properties.
    """
    ver = 'v210'
    bgps = read_cat('bgps_v210_evo').set_index(ver + 'cnum')
    dpdf_props = read_cat('bgps_v210_dpdf_props').set_index(ver + 'cnum')
    omni = read_dpdf(v=2)
    xdist = np.linspace(omni[2].data[0][2],
                        omni[2].data[0][1] * omni[2].data[0][0],
                        omni[2].data[0][0])

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


class PropCollection(object):
    """
    Collection of all DPDFs and their properties. The class also contains
    methods for calculating Monte Carlo sampled properties.
    """
    ver = 'v210'
    bgps = read_cat('bgps_v210_evo').set_index(ver + 'cnum')
    dpdf_props = read_cat('bgps_v210_dpdf_props').set_index(ver + 'cnum')
    omni = read_dpdf(v=2)
    all_props = []

    def __init__(self, posteriors, ktemps):
        self.posteriors = posteriors
        self.ktemps = ktemps
        # get stage

    def construct_dpdfs(self):
        for cnum, post in self.posteriors.iteritems():
            # assign stage
            # append evo props
            # append stage group
            #self.all_props.append(clump_prop)
            pass


class TempDistribs(object):
    temp_cols = ['nh3_gbt_tkin', 'mol_nh3_tk']
    err_cols = ['nh3_gbt_tkin_err', 'mol_nh3_tk_err']

    def __init__(self, stages):
        self.stages = stages

    def get_good_tkins(self):
        for stage in self.stages:
            # NH3 Obs from Rosolowsky
            gbt_mask = (stage['nh3_gbt_snr11'] > 5) & \
                       (stage['nh3_gbt_snr22'] > 5)
            gbt = stage[gbt_mask][['nh3_gbt_tkin', 'nh3_gbt_tkin_err']]
            # NH3 Obs collected from Shirley et al. (2013)
            mol_mask = stage['mol_nh3_tk'].notnull()
            mol = stage[mol_mask][['mol_nh3_tk', 'mol_nh3_tk_err']]
            temps = gbt.merge(mol, how='outer')
            # Average temperatures
            temps['tk'] = temps[self.temp_cols].min(axis=1)
            temps['tk_err'] = temps[self.err_cols].min(axis=1)
            both_mask = gbt_mask & mol_mask
            # Hack for weighted mean
            temps['tk'][both_mask] = \
                (temps['nh3_gbt_tkin'] / temps['nh3_gbt_tkin_err'] +
                 temps['mol_nh3_tk'] / temps['mol_nh3_tk_err']) / \
                (1 / temps['nh3_gbt_tkin_err'] + 1 / temps['mol_nh3_tk_err'])
            temps['tk_err'][both_mask] = \
                np.sqrt(temps['nh3_gbt_tkin_err']**2 +
                        temps['mol_nh3_tk_err']**2) / 2.
            self.good_tkins = temps

    def to_distrib(self):
        self.tk_fn = []
        for stage in self.stages:
            stage_temps = self.good_tkins.ix[stage.index, ['tk', 'tk_err']]
            tk_samples = np.random.normal(stage_temps['tk'].values,
                                          stage_temps['tk_err'].values, 100)
            tk_samples.sort()
            kernel = gaussian_kde(tk_samples)
            kernel.set_bandwidth(stage_temps['tk_err'].median())
            self.tk_fn.append(kernel)
