"""
=========================
DPDF Monte Carlo Sampling
=========================

Provides methods for Monte Carlo sampling DPDF collections.  Also calculate
composite, posterior DPDFs for nodes in PPV-groups.

"""
from __future__ import division
import math
import numpy as np
import pandas as pd
from besl import (catalog, mathf, dpdf_calc)
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from scipy.special import erf


class Prop(object):
    """
    Container object for a BGPS source property.
    """
    ver = 'v210'
    bgps = catalog.read_cat('bgps_v210_evo').set_index(ver + 'cnum')
    dpdf_props = catalog.read_cat('bgps_v210_dpdf_props').set_index(ver + 'cnum')
    omni = catalog.read_dpdf(v=2)
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
        self.y[self.y == 0] = 1e-16  # can't have zero error
        self.nsamples = nsamples


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


class SingleNormalSampler(SamplerBase):
    def draw(self):
        return np.random.normal(self.x, self.y, self.nsamples)


class DistSampler(SamplerBase):
    def __init__(self, data, nsamples):
        super(DistSampler, self).__init__(data, nsamples)
        self.dx = self.x[1] - self.x[0]
        y_cum = np.cumsum(self.y) / np.sum(self.y)
        # Max of `y_cum` is 1, so only out of bounds values should be less
        # than the mininum
        self.y_int = interp1d(y_cum, self.x, bounds_error=False,
            fill_value=y_cum.min())

    def draw(self):
        """
        Take a random draw from the interpolated distribution.

        Returns
        -------
        samples : np.array
        """
        return self.y_int(np.random.uniform(size=self.nsamples)) + self.dx / 2.


class TempSampler(object):
    def __init__(self, cat, nsamples):
        self.cat = cat[cat.nh3_tkin.notnull()]
        self.stages, _ = dpdf_calc.evo_stages(self.cat)
        self.ns = len(self.stages)
        self.nsamples = nsamples
        self.samplers = [NormalSampler((df['nh3_tkin'].values,
                                        df['nh3_tkin_err'].values), nsamples)
                         for df in self.stages]
        self.samples = [s.draw().flatten() for s in self.samplers]
        self.samples = [s[(s > 0) & (s < 200)] for s in self.samples]

    def draw(self, ii):
        return np.random.choice(self.samples[ii], size=self.nsamples)


class DepSampler(object):
    distx = np.arange(1000, dtype=float) * 20. + 20.

    def __init__(self, cat, nsamples):
        assert cat.index.name == 'v210cnum'
        self.cat = cat
        self.nsamples = nsamples
        print ':: Read in data'
        cat = cat.copy()
        self.posts = catalog.read_pickle('ppv_dpdf_posteriors')
        self.dix = {k: v for k, v in self.posts.items()
                    if k in cat.query('10 < glon_peak < 65').index}
        self.stages, _ = dpdf_calc.evo_stages(cat)
        self.ns = len(self.stages)
        self.stage_ix = [df.index for df in self.stages]


class SurfSampler(DepSampler):
    to_solar = 4.788392e3  # (g cm^-2)^-1 Msun pc^-2

    def __init__(self, cat, nsamples, size='full'):
        """
        size : str, Default 'full'
            Size and flux density definition to use for mass surface density
            calculator.
                'full' -> Full Bolocat labelmask equivalent size
                'fwhm' -> FWHM equivalent size
                'peak' -> Peak solid angle (5-pix)
        """
        super(SurfSampler, self).__init__(cat, nsamples)
        self.size = size
        print ':: Sampling fluxes'
        self.fluxes = self.get_fluxes(size)
        self.sangles = self.get_sangles(size)
        print ':: Sampling temperatures (normal)'
        self.tkins = NormalSampler((cat.nh3_tkin.values,
                                    cat.nh3_tkin_err.values), nsamples).draw()
        self.good_tk = cat[cat.nh3_tkin.notnull()].index
        print ':: Sampling temperatures (stages)'
        self.tkin_sampler = TempSampler(cat, nsamples)

    def get_fluxes(self, size):
        if size == 'full':
            return NormalSampler((self.cat['flux'].values,
                                  self.cat['err_flux'].values),
                                 self.nsamples).draw()
        elif size == 'fwhm':
            return NormalSampler((self.cat['fwhm_flux'].values,
                                  self.cat['err_fwhm_flux'].values),
                                 self.nsamples).draw()
        elif size == 'peak':
            return NormalSampler((self.cat['peak_flux'].values,
                                  self.cat['err_peak_flux'].values),
                                 self.nsamples).draw()
        else:
            raise ValueError('Invalid size definition: {0}'.format(size))

    def get_sangles(self, size):
        if size == 'full':
            return self.cat['sangle'].values
        elif size == 'fwhm':
            return self.cat['fwhm_sangle'].values
        elif size == 'peak':
            self.cat['peak_sangle'] = 7.2**2
            return self.cat['peak_sangle'].values
        else:
            raise ValueError('Invalid size definition: {0}'.format(size))

    def draw(self):
        surf = np.empty((self.cat.index.shape[0], self.nsamples), dtype=float)
        surf[:,:] = np.nan
        ovix = self.cat.query('10 < glon_peak < 65').index
        nix = len(ovix)
        for ii, cix in enumerate(ovix):
            print '-- ', ii + 1, ' / ', nix
            csangle = self.sangles[cix - 1]
            cflux = self.fluxes[cix - 1]
            if cix in self.good_tk:
                ctkin = self.tkins[cix - 1]
            else:
                for jj, stix in enumerate(reversed(self.stage_ix)):
                    if cix in stix:
                        break
                else:
                    raise ValueError('Not in stage, cnum: {0}'.format(cix))
                ctkin = self.tkin_sampler.draw(self.ns - 1 - jj)
            surf[cix - 1] = self.calc_surf(cflux, csangle, ctkin)
        return surf

    @staticmethod
    def calc_surf(flux, csangle, tkin):
        return 39.25 * flux * np.pi / csangle * (np.exp(13.08 / tkin) - 1)


class RadiusSampler(DepSampler):
    def __init__(self, cat, nsamples, use_fwhm=False):
        super(RadiusSampler, self).__init__(cat, nsamples)
        self.use_fwhm = use_fwhm
        if use_fwhm:
            self.angle = cat['fwhm_eqangled'].values
        else:
            self.angle = cat['eqangled'].values

    def draw(self):
        dix = self.dix
        radii = np.empty((self.cat.shape[0], self.nsamples), dtype=float)
        radii[:,:] = np.nan
        for ii, cix in enumerate(dix):
            print '-- ', ii + 1, ' / ', len(dix)
            cdist = DistSampler((self.distx, self.posts[cix]), self.nsamples).draw()
            cangle = self.angle[cix - 1]
            radii[cix - 1] = self.calc_radius(cangle, cdist)
        return radii

    @staticmethod
    def calc_radius(angle, dist):
        # `angle` must be in arcsec
        # `dist` must be in pc
        arc2rad = 4.848137e-6
        return dist * arc2rad * angle


class CenDensSampler(DepSampler):
    def __init__(self, cat, nsamples):
        super(CenDensSampler, self).__init__(cat, nsamples)
        self.rs = RadiusSampler(cat, nsamples, use_fwhm=True)
        self.ss = SurfSampler(cat, nsamples, size='peak')

    def draw(self):
        radii = self.rs.draw()
        surfs = self.ss.draw()
        teq = np.array([self.cat['eqangled'].values]).T
        feq = np.array([self.cat['fwhm_eqangled'].values]).T
        mu = 2.8  # ref Dunham et al. 2011
        a1 = np.sqrt(np.log(2) / np.pi) / (mu * 1.672622e-24)  # g -> [H2]
        a2 = 2 * 3.0856776e18  # pc -> cm
        a3 = np.sqrt(np.log(2))
        return a1 * surfs / (a2 * radii) / erf(a3 * teq / feq)


class AvgDensSampler(DepSampler):
    def __init__(self, cat, nsamples, use_fwhm=False):
        super(AvgDensSampler, self).__init__(cat, nsamples)
        self.rs = RadiusSampler(cat, nsamples, use_fwhm)
        self.ms = MassSampler(cat, nsamples, use_fwhm)

    def draw(self):
        radii = self.rs.draw()
        masses = self.ms.draw()
        mu = 2.8  # ref Dunham et al. 2011
        mp = 1.672622e-24
        a1 = 6.7702543e-23 / (mu * mp)  # Msun / pc^3 -> [H2] / cm^3
        a2 = 4 * np.pi / 3.
        return a1 * masses / (a2 * radii**3)


class CylDensSampler(DepSampler):
    def __init__(self, cat, nsamples):
        super(CylDensSampler, self).__init__(cat, nsamples)
        self.rst = RadiusSampler(cat, nsamples, use_fwhm=False)
        self.rsf = RadiusSampler(cat, nsamples, use_fwhm=True)
        self.ms = MassSampler(cat, nsamples, use_fwhm=True)

    def draw(self):
        tot_radii = self.rst.draw()
        fwhm_radii = self.rsf.draw()
        masses = self.ms.draw()
        mu = 2.8  # ref Dunham et al. 2011
        mp = 1.672622e-24
        a1 = 6.7702543e-23 / (mu * mp)  # Msun / pc^3 -> [H2] / cm^3
        return a1 * masses / (2 * np.pi * fwhm_radii**2 * tot_radii)


class VirSampler(DepSampler):
    pix = 1.8
    pix_err = 0.4

    def __init__(self, cat, nsamples):
        cat = self.clean_cat(cat)
        super(VirSampler, self).__init__(cat, nsamples)
        self.vs = NormalSampler((cat['nh3_gbt_sigmav'].values,
                                 cat['nh3_gbt_sigmav_err'].values), nsamples)
        self.rs = RadiusSampler(cat, nsamples, use_fwhm=False)
        self.ms = MassSampler(cat, nsamples, use_fwhm=False)
        nn = self.cat.shape[0]
        self.ps = NormalSampler((self.pix * np.ones(nn),
                                 self.pix_err * np.ones(nn)), nsamples)

    def clean_cat(self, cat):
        bad_velos = cat.query('nh3_gbt_sigmav_err == 0').index
        cat.loc[bad_velos, 'nh3_gbt_sigmav'] = np.nan
        cat.loc[bad_velos, 'nh3_gbt_sigmav_err'] = np.nan
        return cat

    def draw(self):
        velos = self.vs.draw()
        radii = self.rs.draw()
        masses = self.ms.draw()
        p = self.ps.draw()
        a1 = (1 - p / 3.) / (1 - 2 * p / 5.)
        a1[(a1 <= 0) | (a1 > 5)] = np.nan
        bigG = 0.0043021135  # km^2 s^-2 pc Msun^-1
        virs = (5 / (8 * np.log(2))) * (2.35482 * velos)**2 * radii / (a1 * bigG * masses)
        # clean virs
        virs[(virs < 1e-3) | (1e3 < virs)] = np.nan
        return virs


class MassSampler(DepSampler):
    def __init__(self, cat, nsamples, use_fwhm=False):
        super(MassSampler, self).__init__(cat, nsamples)
        self.use_fwhm = use_fwhm
        # flux samples, index offset of -1
        print ':: Sampling fluxes'
        self.fluxes = self.get_fluxes()
        print ':: Sampling temperatures (normal)'
        self.tkins = NormalSampler((cat.nh3_tkin.values,
                                    cat.nh3_tkin_err.values), nsamples).draw()
        self.good_tk = cat[cat.nh3_tkin.notnull()].index
        print ':: Sampling temperatures (stages)'
        self.tkin_sampler = TempSampler(cat, nsamples)

    def get_fluxes(self):
        if self.use_fwhm:
            return NormalSampler((self.cat['fwhm_flux'].values,
                                  self.cat['err_fwhm_flux'].values),
                                 self.nsamples).draw()
        else:
            return NormalSampler((self.cat['flux'].values,
                                  self.cat['err_flux'].values),
                                 self.nsamples).draw()

    def draw(self):
        dix = self.dix
        masses = np.empty((self.cat.shape[0], self.nsamples), dtype=float)
        masses[:,:] = np.nan
        for ii, cix in enumerate(dix):
            print '-- ', ii + 1, ' / ', len(dix)
            cdist = DistSampler((self.distx, self.posts[cix]), self.nsamples).draw()
            cflux = self.fluxes[cix - 1]
            if cix in self.good_tk:
                ctkin = self.tkins[cix - 1]
            else:
                for jj, stix in enumerate(reversed(self.stage_ix)):
                    if cix in stix:
                        break
                else:
                    raise ValueError('Not in stage, cnum: {0}'.format(cix))
                ctkin = self.tkin_sampler.draw(self.ns - 1 - jj)
            masses[cix - 1] = self.calc_mass(ctkin, cflux, cdist)
        return masses

    @staticmethod
    def calc_mass(tkin, flux, dist):
        # dist must be in pc
        return 14.067 * (np.exp(13.08 / tkin) - 1) * flux * (dist * 1e-3)**2


class FreeFallSampler(DepSampler):
    def __init__(self, cat, nsamples, use_fwhm=True):
        super(FreeFallSampler, self).__init__(cat, nsamples)
        self.use_fwhm = use_fwhm
        self.ms = MassSampler(self.cat, nsamples, use_fwhm=True)

    def get_volumes(self):
        dix = self.dix
        vols = np.empty((self.cat.shape[0], self.nsamples), dtype=float)
        vols[:,:] = np.nan
        for ii, cix in enumerate(dix):
            cdist = DistSampler((self.distx, self.posts[cix]), self.nsamples).draw()
            vols[cix - 1] = self.calc_volume(cdist, self.ms.cat.loc[cix, 'fwhm_eqangled'])
        return vols

    def draw(self):
        print ':: Sampling masses'
        self.masses = self.ms.draw()
        print ':: Sampling volumes'
        self.volumes = self.get_volumes()
        print ':: Calculating freefall times'
        return self.calc_freefall(self.masses, self.volumes)

    @staticmethod
    def calc_volume(dist, angle):
        # `dist` in parsec
        # `angle` in arcsec
        arcsec2rad = 4.8481368e-6
        return 4 * np.pi / 3. * (arcsec2rad * angle * dist / 2.)**3

    @staticmethod
    def calc_freefall(mass, volume):
        msun_per_pc3 = 6.7702543e-20  # in kg / m^3
        return 0.002680274 / np.sqrt(msun_per_pc3 * mass / volume)  # in yr


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
        cat['robit_yso'] = cat.robit_f.copy().replace(np.nan, 0).astype(int)
        cat['robit_agb'] = cat.eval('robit_n - robit_f').replace(np.nan, 0).astype(int)
        cat['robit_yso_r'] = 0
        cat['robit_agb_r'] = 0
        cat['ir_f'] = 0
        self.cat = cat

    def draw_stages(self):
        cat = self.cat  # still the same dataframe in memory
        cat.robit_yso_r = 0
        cat.robit_agb_r = 0
        cat.ir_f = 0
        y2y = np.random.binomial(cat.robit_yso, p=self.cfrac_yso)
        a2a = np.random.binomial(cat.robit_agb, p=self.cfrac_agb)
        # compliments are the resampled to the other category
        cat.robit_yso_r = cat.eval('@y2y + robit_agb - @a2a')  # y2y + a2y
        cat.robit_agb_r = cat.eval('@a2a + robit_yso - @y2y')  # a2a + y2a
        # Evo stages will select based on `ir_f` and `sf_f`.
        # Because agb's are not a positive indicator, only the resampled YSO's
        # count towards the flags.
        cat.loc[(cat.robit_yso_r > 0) |
                (cat.red_msx_f > 0) |
                (cat.ego_n > 0), 'ir_f'] = 1
        cat.loc[(cat.ir_f == 1) |
                (cat.h2o_f == 1) |
                (cat.ch3oh_f == 1) |
                (cat.uchii_f == 1), 'sf_f'] = 1
        stages, labels = dpdf_calc.evo_stages(bgps=cat)
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
    bgps = catalog.read_cat('bgps_v210_evo').set_index(ver + 'cnum')
    dpdf_props = catalog.read_cat('bgps_v210_dpdf_props').set_index(ver + 'cnum')
    omni = catalog.read_dpdf(v=2)
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
                temps = mathf.pd_weighted_mean(temps, self.temp_cols,
                                               self.err_cols,
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


