"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""

import numpy as np
from .catalog import read_cat, read_dpdf
from .mathf import pd_weighted_mean
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


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
        # check whether has dpdf
        # pdf
        # kinetic temperature
        # initiate sampler
        pass


class Sampler(object):
    """
    Simple Monte Carlo sampler.
    """
    def __init__(self, data, normal=False):
        """
        Initialize the sampler.

        Parameters
        ----------
        data : tuple
            Data can be a length 2 tuple containing either:
                (mu, sigma) : (number, number)
                    To draw random normal deviates.
                (x, y) : (array, array)
                    x & y values of a distribution to interpolate
        normal : bool
            Whether to use normal distribution
        """
        assert len(data) == 2
        self.data = data
        self.normal = normal

    def draw(self, nsamples=1e3):
        """
        Take a random draw from the distribution.

        Parameters
        ----------
        nsamples : number
            Number of samples to draw

        Returns
        -------
        samples : np.array
        """
        # Normal deviates
        if self.normal:
            return np.random.normal(self.data[0], self.data[1], nsamples)
        # Random draw from general distribution
        else:
            x, y = self.data
            dx = x[1] - x[0]
            y_cum = np.cumsum(y) / np.sum(y)
            # Max of `y_cum` is 1, so only out of bounds values should be less
            # than the mininum
            y_int = interp1d(y_cum, x, bounds_error=False,
                             fill_value=y_cum.min())
            samples = y_int(np.random.uniform(size=nsamples))
            samples += dx / 2.
            return samples


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
    v = 'v210'
    temp_cols = ['nh3_gbt_tkin', 'mol_nh3_tk']
    err_cols = ['nh3_gbt_tkin_err', 'mol_nh3_tk_err']
    out_cols = ['tk', 'tk_err']

    def __init__(self, stages):
        self.stages = stages
        self.good_tkins = []
        self.tk_fn = []
        self.get_good_tkins()
        self.to_distrib()

    def get_good_tkins(self, weight=True):
        for stage in self.stages:
            # NH3 Obs from Rosolowsky
            gbt_mask = (stage['nh3_gbt_snr11'] > 5) & \
                       (stage['nh3_gbt_snr22'] > 3)
            gbt = stage[gbt_mask][['v210cnum', 'nh3_gbt_tkin',
                                   'nh3_gbt_tkin_err']]
            # NH3 Obs collected from Shirley et al. (2013)
            mol_mask = stage['mol_nh3_tk'].notnull()
            mol = stage[mol_mask][[self.v + 'cnum', 'mol_nh3_tk',
                                   'mol_nh3_tk_err']]
            temps = gbt.merge(mol, how='outer', on=self.v + 'cnum')
            # Weighted mean temperatures
            if weight:
                temps = pd_weighted_mean(temps, self.temp_cols, self.err_cols,
                                         out_cols=self.out_cols)
            # Use GBT over collected
            else:
                temps['tk'] = temps[self.temp_cols[1]]
                temps['tk'][temps[self.temp_cols[0]].notnull()] = \
                    temps[self.temp_cols[0]]
                temps['tk_err'] = temps[self.err_cols[1]]
                temps['tk_err'][temps[self.err_cols[0]].notnull()] = \
                    temps[self.temp_cols[0]]
            self.good_tkins.append(temps)

    def to_distrib(self):
        for temps in self.good_tkins:
            kernel = gaussian_kde(temps[temps['tk'].notnull()]['tk'])
            kernel.set_bandwidth(temps['tk_err'].median())
            self.tk_fn.append(kernel)


