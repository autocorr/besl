"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""
# TODO add in PPV grouping properties

import numpy as np
import pandas as pd
from besl.catalog import (read_cat, read_dpdf)
from besl.mathf import (pd_weighted_mean)
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d


class Prop(object):
    """
    Container object for a BGPS source property.
    """
    ver = 'v210'
    bgps = read_cat('bgps_v210_evo').set_index(ver + 'cnum')
    dpdf_props = read_cat('bgps_v210_dpdf_props').set_index(ver + 'cnum')
    omni = read_dpdf(v=2)
    xdist = np.linspace(omni[2].data[0][2],
                        omni[2].data[0][1] * omni[2].data[0][0],
                        omni[2].data[0][0])

    def __init__(self, data, normal=None):
        self.data = data
        self.normal = normal
        # check whether has dpdf
        # pdf
        # kinetic temperature
        # initiate sampler
        self.sampler = Sampler(self.data, normal=self.normal)
        pass

    def draw(self):
        pass


class SamplerBase(object):
    """
    Simple Monte Carlo sampler base class.
    """
    def __init__(self, data, nsamples):
        """
        Initialize the sampler.

        Parameters
        ----------
        data : tuple, (array-like, array-like)
        nsamples : number
            Number of samples to draw
        """
        assert len(data) == 2
        self.x = data[0]
        self.y = data[1]
        self.nsamples = nsamples
        self.shape = (len(self.x), nsamples)


class NormalSampler(SamplerBase):
    def draw(self):
        """
        Take a random draw from a Gaussian distribution.

        Returns
        -------
        samples : np.array
        """
        return np.array([
            np.random.normal(ii, jj, self.nsamples)
                for ii, jj in
            np.array([self.x, self.y]).T])


class DistSampler(SamplerBase):
    def draw(self):
        """
        Take a random draw from the interpolated distribution.

        Returns
        -------
        samples : np.array
        """
        dx = self.x[1] - self.x[0]
        y_cum = np.cumsum(self.y) / np.sum(self.y)
        # Max of `y_cum` is 1, so only out of bounds values should be less
        # than the mininum
        y_int = interp1d(y_cum, self.x, bounds_error=False,
            fill_value=y_cum.min())
        return y_int(np.random.uniform(size=self.nsamples)) + dx / 2.


class BoolResampler(object):
    def __init__(self, ix1, ix2, cfrac1=0.5, cfrac2=0.5):
        """
        Resample two lists of indices with what fractions to move to the other
        group.

        Parameters
        ----------
        ix1 : pandas.core.index
        ix2 : pandas.core.index
        cfrac1 : number, default 0.5
            Fraction of ix1 to move into ix2
        cfrac2 : number, default 0.5
            Fraction of ix2 to move into ix1
        """
        self.ix1 = ix1
        self.ix2 = ix2
        self.cfrac1 = cfrac1
        self.cfrac2 = cfrac2

    def resample(self):
        contam1 = self._draw(self.ix1, self.cfrac1)
        contam2 = self._draw(self.ix2, self.cfrac2)
        rix1 = self.ix1.drop(contam1).append(contam2)
        rix2 = self.ix2.drop(contam2).append(contam1)
        return rix1, rix2

    @staticmethod
    def _draw(ix, cfrac):
        contam = np.random.choice(ix, size=round(len(ix) * cfrac), replace=False)
        return pd.Index(contam)


class IrResampler(object):
    cfrac_yso = 0.5
    cfrac_agb = 0.6

    def __init__(self, cat=None):
        """
        Calculate a resampling of the evo stages for the IR Robitaille
        indicator.  Uses the binomial distribution over the number of YSO or AGB
        counts within a clump given a "contamination fraction".

        Parameters
        ----------
        cat : pandas.DataFrame

        Attributes
        ----------
        cfrac_yso : float, default 0.5
            Probability of a YSO being a YSO
        cfrac_agb : float, default 0.6
            Probability of a AGB being an AGB
        """
        if cat is None:
            cat = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
        else:
            assert cat.index.name == 'v210cnum'
        cat['robit_yso'] = cat.robit_f.copy()
        cat['robit_agb'] = cat.eval('robit_n - robit_f')
        cat['robit_yso_r'] = 0
        cat['robit_agb_r'] = 0
        cat['robit_agb_r'] = 0
        self.cat = cat

    def resample(self):
        self.cat.robit_yso_r = 0
        self.cat.robit_agb_r = 0
        new_yso = self.cat.robit_yso.apply(np.random.binomial,
                                              p=self.cfrac_yso)
        new_agb = self.cat.robit_yso - yso_to_agb
        y2y = self.cat.robit_yso.apply(np.random.binomial, p=self.cfrac_yso)
        y2a = self.cat.robit_yso - y2y  # compliment are the resampled agb's
        a2a = self.cat.robit_agb.apply(np.random.binomial, p=self.cfrac_agb)
        a2y = self.cat.robit_agb - a2a
        self.cat.robit_yso_r = y2y + a2y
        self.cat.robit_agb_r = a2a + y2a

    def draw_stages(self):
        self.resample()
        self.cat.ir_f = 0
        # Evo stages will select based on `ir_f` and `sf_f`.
        # Because agb's are not a positive indicator, only the resampled YSO's
        # count towards the flags.
        self.cat.loc[(rcat.robit_yso_r > 0) |
                     (rcat.red_msx_f > 0) |
                     (rcat.ego_n > 0), 'ir_f'] = 1
        self.cat.loc[(rcat.ir_f == 1) |
                     (rcat.h2o_f == 1) |
                     (rcat.ch3oh_f == 1) |
                     (rcat.uchii_f == 1), 'sf_f'] = 1
        stages, labels = dpdf_calc.evo_stages(bgps=self.cat)
        return stages, labels


class MedianStats(object):
    cols = ['med1', 'med2', 'ks_mu', 'ks_p']
    normal = True

    def __init__(self, cat, prop, eprop, Sampler, nsamples):
        """
        Calculate median statistics from a catalog for a property and it's
        uncertainty.

        Parameters
        ----------
        cat : pd.DataFrame
        prop : string or array-like
            Property to MC
        eprop : string or array-like
            Uncertainty of property. Passed to `Sampler`.
        Sampler : Sampler
            Sampler class to draw random samples from `prop` `eprop`.
        nsamples : number
            Number of samples to draw.
        """
        self.cat = cat
        self.prop = prop
        self.eprop = eprop
        self.Sampler = Sampler
        self.nsamples = nsamples

    def draw(self):
        x, y = self.cat.loc[ii, [self.prop, self.eprop]]
        sampler = self.Sampler((prop, eprop), self.nsamples)
        return pd.DataFrame(data=sampler.draw(), index=self.cat.index,
                               columns=np.arange(self.nsamples))

    def diff_stats(ix1, ix2, samples):
        stats = pd.DataFrame(index=np.arange(samples.shape[1]),
                             columns=self.cols)
        for col in samples.columns:
            s1 = samples.loc[ix1, col]
            s2 = samples.loc[ix2, col]
            med1 = s1.median()
            med2 = s2.median()
            ks_mu, ks_p = ks_2samp(s1, s2)
            stats.loc[col] = [med1, med2, ks_mu, ks_p]
        stats['meddiff'] = stats.med2 - stats.med1
        return stats


class PropCollection(object):
    """
    Collection of all like properties. The class also contains
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

    def sample(self):
        pass

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
            gbt_mask = (stage['nh3_gbt_snr11'] > 3) & \
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


###
# Specific Observed Quantites
###

# kinetic temperature
# molecular velocity width
# molecular integrated intensity
# flux
# angular diameter
# angular size
# water maser properties

###
# Specific Derived Quantities
###

# distance
# surface mass density
# dust mass
# physical equivalent radius
# physical area
# mass density
# virial parameter
# water maser isotropic luminosity


