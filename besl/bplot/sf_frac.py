#!/usr/bin/env python
# encoding: utf-8
"""
=======================
Star-Formation Fraction
=======================

Plot the starless fraction or protostellar fraction as a function of a property.

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from besl import catalog
from besl import dpdf_calc
from besl import dpdf_mc
from besl import util


plt.rc('font', **{'size':10})


def get_data(clip_lon=True):
    """
    Get combined dataframe of evolutionary and fwhm catalogs.

    Parameters
    ----------
    clip_lon : bool, Default True
        Clip the dataframe to 10 < l < 65 degrees

    Returns
    -------
    df : pd.DataFrame
    """
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    fwhm = catalog.read_cat('bgps_v210_fwhm').set_index('v210cnum')
    fwhm = fwhm.loc[:, 'npix':]
    df = evo.merge(fwhm, left_index=True, right_index=True)
    if clip_lon:
        df = df.query('10 < glon_peak < 65')
    return df


class ObsData(object):
    winr = 25

    def __init__(self, col, xlabel, with_dist=False, bool_proto=False):
        """
        Star formation fraction for observable quantities without MC simulation.

        Parameters
        ----------
        col : str
            Dataframe column name of property
        xlabel : str
            LaTeX string for xlabel of plot
        with_dist : bool, Default False
            Calculate values for the Distance Sample sub-set
        bool_proto : bool, Default False
            Use the boolean flags for protostellar activity

        Attributes
        ----------
        winr : number, Default 25
            Window radius used in boxcar fraction calculation
        """
        # params
        self.col = col
        self.xlabel = xlabel
        self.with_dist = with_dist
        self.bool_proto = bool_proto
        # labels
        self.sf_fcol = 'sf_frac_' + col
        self.sf_tcol = 'sf_tot_' + col
        self.sf_wcol = 'sf_win_' + col
        if bool_proto:
            self.proto_col = 'is_proto'
        else:
            self.proto_col = 'proto_prob'
        # data
        self.df = self.get_df(with_dist=with_dist)
        self.mean_frac = self.df[self.proto_col].mean()
        self.cumix = self.df.index
        self.bins = self.df[col].values
        self.nbins = self.df.shape[0]
        self.xmin, self.xmax = self.bins[0], self.bins[-1]
        # compute
        print ':: Computing {0}'.format(col)
        self.fvals, self.tvals, self.wvals = self.cumfrac(self.df, self.cumix)

    def get_df(self, with_dist=False):
        df = get_data()
        if with_dist:
            posts = catalog.read_pickle('ppv_dpdf_posteriors')
            dix = {k: v for k, v in posts.items()
                   if k in df.query('10 < glon_peak < 65').index}
            df = df.loc[dix.keys()]
        df = df[df[self.col].notnull() & (df[self.col] > 0)].sort(self.col)
        stages, labels = dpdf_calc.evo_stages(bgps=df)
        df['is_proto'] = False
        df.loc[stages[1].index, 'is_proto'] = True
        for col in (self.sf_fcol, self.sf_tcol, self.sf_wcol):
            df[col] = np.nan
        return df

    def cumfrac(self, df, cumix):
        nsamp = df.shape[0]
        for nn, ii in enumerate(cumix):
            nproto = df.loc[:ii, self.proto_col].sum()
            df.loc[ii, self.sf_fcol] = nproto / nn
            df.loc[ii, self.sf_tcol] = nproto / nsamp
            df.loc[ii, self.sf_wcol] = self.window(nn, df, cumix)
        return (df[self.sf_fcol].values,
                df[self.sf_tcol].values,
                df[self.sf_wcol].values)

    def window(self, nn, df, cumix):
        nsamp = df.shape[0]
        # repeats last value off edge
        if nn < self.winr:
            win = cumix.values[:self.winr+1]
        elif nn > nsamp - self.winr:
            win = cumix.values[-self.winr:]
        else:
            win = cumix.values[nn-self.winr:nn+self.winr+1]
        return df.loc[win, self.proto_col].mean()


class DistData(object):
    col = 'dist'
    xlabel = r'${\rm Heliocentric \ Distance} \ \ [{\rm kpc}]$'
    distx = np.arange(1e3) * 20. + 20.
    xmin = distx[0]
    xmax = distx[-1]
    proto_col = 'proto_prob'
    up_dist = 12000  # pc

    def __init__(self):
        print ':: Read in posteriors'
        self.df = get_data()
        self.posts = catalog.read_pickle('ppv_dpdf_posteriors')
        self.dix = {k: v for k, v in self.posts.items()
                    if k in self.df.query('10 < glon_peak < 65').index}
        self.mean_frac = self.df.loc[self.dix.keys(), self.proto_col].mean()
        self.sf_frac = self.calc()
        self.loix, self.nup, self.up_frac = self.calc_dlim()
        self.up_label = self.gen_label()

    def calc(self):
        print ':: Compute protostellar fraction'
        sf_frac = np.zeros(self.distx.shape, dtype=float)
        dpdf_sum = np.zeros(self.distx.shape, dtype=float)
        for ii, dd in self.dix.items():
            sf_frac += self.df.loc[ii, self.proto_col] * dd
            dpdf_sum += dd
        sf_frac /= dpdf_sum
        return sf_frac

    def calc_dlim(self):
        loix = np.argwhere(self.distx == self.up_dist - 20)[0,0]
        nup = sum([1 for k, v in self.dix.items()
                   if self.distx[v.argmax()] > self.up_dist])
        up_frac = nup / len(self.dix)
        return loix, nup, up_frac

    def gen_label(self):
        return (r'$N(d_{\rm ML} >' +
                '{0:2.0f}'.format(self.up_dist / 1e3) +
                r'\ {\rm kpc})=' +
                '{0:1.1f}\% \\ '.format(self.up_frac * 1e2) +
                '({0})'.format(self.nup) +
                '$')


class MassData(object):
    col = 'mass'
    xlabel = r'$M^{\rm Total}_{\rm H_2} \ \ [{\rm M_\odot}]$'
    use_fwhm = False
    nsamples = 1e4
    bins = np.logspace(1, 5, 300)
    xmin = bins[0]
    xmax = bins[-1]
    proto_col = 'proto_prob'

    def __init__(self):
        self.df = get_data(clip_lon=False)
        if self.use_fwhm:
            self.xlabel = r'$M^{\rm FWHM}_{\rm H_2} \ \ [{\rm M_\odot}]$'
            self.col += '_fwhm'
            bad_fwhm = self.df.query('fwhm_flux == 0 | err_fwhm_flux == 0').index
            self.df.loc[bad_fwhm, 'fwhm_flux'] = np.nan
            self.df.loc[bad_fwhm, 'err_fwhm_flux'] = np.nan
        self.ms = dpdf_mc.MassSampler(self.df, self.nsamples, use_fwhm=self.use_fwhm)
        self.mean_frac = self.df.loc[self.ms.dix.keys(), self.proto_col].mean()
        self.sf_frac = self.calc()

    def calc(self):
        masses = self.ms.draw()
        hist = np.zeros(len(self.bins)-1, dtype=float)
        mass_sum = np.zeros(len(self.bins)-1, dtype=float)
        print ':: Compute histogram'
        for ii, dd in self.ms.dix.items():
            mm = masses[ii-1]
            pp = self.df.loc[ii, self.proto_col]
            hh, _ = np.histogram(mm, bins=self.bins)
            hist += hh * pp / self.nsamples
            mass_sum += hh / self.nsamples
        hist /= mass_sum
        return hist


class MassFracData(MassData):
    corr_fact = 1015. / 625.  # DPDF sample correction factor
    mean_frac = 1 / corr_fact

    def calc(self):
        masses = self.ms.draw()
        scc_hist = np.zeros(len(self.bins)-1, dtype=float)
        pro_hist = scc_hist.copy()
        print ':: Compute histogram'
        for ii, dd in self.ms.dix.items():
            mm = masses[ii-1]
            pp = self.df.loc[ii, self.proto_col]
            hh, _ = np.histogram(mm, bins=self.bins)
            scc_hist += hh * (1 - pp) / self.nsamples
            pro_hist += hh * pp / self.nsamples
        return self.corr_fact * scc_hist / pro_hist


class LifeTimeData(MassData):
    maser_col = 'ch3oh_f'

    def calc(self):
        masses = self.ms.draw()
        scc_hist = np.zeros(len(self.bins)-1, dtype=float)
        mas_hist = np.zeros(len(self.bins)-1, dtype=float)
        print ':: Compute histogram'
        for ii, dd in self.ms.dix.items():
            mm = masses[ii-1]
            if self.df.loc[ii, self.maser_col] == 1:
                hh, _ = np.histogram(mm, bins=self.bins)
                mas_hist += hh / self.nsamples
            elif self.df.loc[ii, self.proto_col] <= 0.4:
                hh, _ = np.histogram(mm, bins=self.bins)
                scc_hist += hh / self.nsamples
            else:
                pass
        return scc_hist / mas_hist


class Plot(object):
    fontsize = 9

    def __init__(self, od):
        self.od = od

    def make_fig(self):
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.set_xscale('log')
        ax.set_xlabel(self.od.xlabel)
        ax.set_xlim(self.od.xmin, self.od.xmax)
        ax.set_ylim(0, 1.1)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25)
        return fig, ax


class ObsPlot(Plot):
    def plot_frac(self):
        fig, ax = self.make_fig()
        ax.set_ylabel(r'$R_{\rm proto}(X \leq x)$')
        ax.hlines(1, self.od.xmin, self.od.xmax, linestyles='dashed',
                  colors='0.5')
        ax.plot(self.od.bins, self.od.fvals, 'k-', drawstyle='steps')
        util.savefig(self.od.sf_fcol, close=True)

    def plot_tot(self):
        fig, ax = self.make_fig()
        ax.set_ylabel(r'${\rm CDF(Protostellar)}$')
        ax.hlines(self.od.tvals[-1], self.od.xmin, self.od.xmax,
                  linestyles='dashed', colors='0.5')
        ax.plot(self.od.bins, self.od.tvals, 'k-', drawstyle='steps')
        util.savefig(self.od.sf_tcol, close=True)

    def plot_win(self, dod=None, window_label=False):
        fig, ax = self.make_fig()
        if window_label:
            ax.annotate(r'${\rm Window} = ' + str(2 * self.od.winr) + '$',
                        xy=(0.725, 0.1), xycoords='axes fraction',
                        fontsize=self.fontsize)
        ax.set_ylabel(r'$R_{\rm proto}$')
        # if over-plot with the distance sample
        if dod is not None:
            ax.annotate('$({0})$'.format(dod.nbins), xy=(0.860, 0.1),
                        xycoords='axes fraction', fontsize=self.fontsize,
                        color='0.5')
            ax.hlines(dod.mean_frac, self.od.xmin, self.od.xmax,
                      linestyles='dotted', colors='0.5')
            ax.plot(dod.bins, dod.wvals, color='0.5', linestyle='solid',
                    drawstyle='steps')
            ax.plot(dod.bins[:dod.winr], dod.wvals[:dod.winr],
                    color='0.75', linestyle='solid', drawstyle='steps')
            ax.plot(dod.bins[-dod.winr:], dod.wvals[-dod.winr:],
                    color='0.75', linestyle='solid', drawstyle='steps')
            ax.vlines(dod.bins, 1.01, 1.040, colors='0.35',
                      linestyles='solid', linewidth=0.1)
        ax.annotate('$N={0}$'.format(self.od.nbins), xy=(0.675, 0.1),
                    xycoords='axes fraction', fontsize=self.fontsize)
        ax.hlines(1, self.od.xmin, self.od.xmax, linestyles='dashed',
                  colors='0.5')
        ax.hlines(self.od.mean_frac, self.od.xmin, self.od.xmax,
                  linestyles='dotted', colors='black')
        ax.vlines(self.od.bins, 1.040, 1.07, colors='black', linestyles='solid',
                  linewidth=0.1)
        ax.plot(self.od.bins, self.od.wvals, 'k-', drawstyle='steps')
        ax.plot(self.od.bins[:self.od.winr], self.od.wvals[:self.od.winr],
                color='0.5', linestyle='solid', drawstyle='steps')
        ax.plot(self.od.bins[-self.od.winr:], self.od.wvals[-self.od.winr:],
                color='0.5', linestyle='solid', drawstyle='steps')
        util.savefig('{0}_{1}'.format(self.od.sf_wcol, self.od.winr),
                     close=True)


class DistPlot(Plot):
    def plot_frac(self, up_label=False):
        fig, ax = self.make_fig()
        if up_label:
            ax.annotate(self.od.up_label, xy=(0.47, 0.1),
                        xycoords='axes fraction', fontsize=self.fontsize)
            ax.hlines(0.165, 7.5, 9, colors='red')
        ax.annotate(r'$N={0}$'.format(len(self.od.dix)),
                     xy=(0.675, 0.1), xycoords='axes fraction',
                    fontsize=self.fontsize)
        ax.hlines(1, self.od.xmin / 1e3, self.od.xmax / 1e3,
                  linestyles='dashed', colors='0.5')
        ax.hlines(self.od.mean_frac, self.od.xmin / 1e3, self.od.xmax / 1e3,
                  linestyles='dotted', colors='0.5')
        ax.plot(self.od.distx / 1e3, self.od.sf_frac, 'k-')
        ax.plot(self.od.distx[self.od.loix:] / 1e3,
                self.od.sf_frac[self.od.loix:], 'r-')
        ax.set_xscale('linear')
        ax.set_ylabel(r'$R_{\rm proto}$')
        ax.set_xlim(self.od.xmin / 1e3, self.od.xmax / 1e3)
        ax.set_xticks(range(0, 21, 1), minor=True)
        util.savefig('sf_frac_' + self.od.col, close=True)


class MassPlot(Plot):
    def plot_frac(self):
        fig, ax = self.make_fig()
        ax.annotate(r'$N={0}$'.format(len(self.od.ms.dix)),
                     xy=(0.675, 0.1), xycoords='axes fraction',
                    fontsize=self.fontsize)
        ax.hlines(1, self.od.xmin, self.od.xmax,
                  linestyles='dashed', colors='0.5')
        ax.hlines(self.od.mean_frac, self.od.xmin, self.od.xmax,
                  linestyles='dotted', colors='0.5')
        ax.plot(self.od.bins[:-1], self.od.sf_frac, 'k-', drawstyle='steps')
        ax.set_ylabel(r'$R_{\rm proto}$')
        ax.set_xlim(self.od.xmin, self.od.xmax)
        util.savefig('sf_frac_' + self.od.col, close=True)


class MassFracPlot(Plot):
    def plot_frac(self):
        fig, ax = self.make_fig()
        ax.set_ylim([0, 2.5])
        ax.annotate(r'$N={0}$'.format(len(self.od.ms.dix)),
                     xy=(0.675, 0.8), xycoords='axes fraction',
                    fontsize=self.fontsize)
        ax.hlines(1, self.od.xmin, self.od.xmax,
                  linestyles='dashed', colors='0.5')
        ax.hlines(self.od.mean_frac, self.od.xmin, self.od.xmax,
                  linestyles='dotted', colors='0.5')
        ax.plot(self.od.bins[:-1], self.od.sf_frac, 'k-', drawstyle='steps')
        ax.set_ylabel(r'$N_{\rm SCC} / N_{\rm proto}$')
        ax.set_xlim(self.od.xmin, self.od.xmax)
        util.savefig('sf_frac_' + self.od.col, close=True)


class LifeTimePlot(Plot):
    hi_life = 4.0e4
    lo_life = 2.5e4
    max_ix = 225

    def plot_frac(self):
        fig, ax = self.make_fig()
        ax.annotate(r'$N={0}$'.format(len(self.od.ms.dix)),
                     xy=(0.675, 0.85), xycoords='axes fraction',
                    fontsize=self.fontsize)
        bins = self.od.bins[:self.max_ix]
        samples = self.od.sf_frac[:self.max_ix]
        ax.fill_between([1e2, 4e2], 1e3, 1e7, facecolor='0.9',
                        edgecolor='none',)
        ax.fill_between(bins, self.lo_life * samples, self.hi_life * samples,
                        facecolor='0.5', edgecolor='none')
        ax.plot(bins, self.lo_life * samples, 'k-', drawstyle='steps')
        ax.plot(bins, self.hi_life * samples, 'k-', drawstyle='steps')
        ax.plot([bins[-1], bins[-1]],
                [self.lo_life * samples[-1], self.hi_life * samples[-1]], 'k-')
        ax.plot([1.2e4, 8e4], [9e3, 9e3], 'k-')
        ax.plot([3.2e4], [7.3e3], 'kv', ms=5)
        ax.set_yscale('log')
        ax.set_ylim([1e3, 1e7])
        ax.set_ylabel(r'$\tau_{\rm SCC} \ \ [{\rm yr}]$')
        ax.set_xlim(1e2, self.od.xmax)
        util.savefig('lifetime_' + self.od.col, close=True)


def plot_flux_window():
    od = ObsData('flux', r'$S^{\rm Total}_{1.1} \ \ [{\rm Jy}]$')
    dod = ObsData('flux', r'$S^{\rm Total}_{1.1} \ \ [{\rm Jy}]$',
                  with_dist=True)
    op = ObsPlot(od)
    op.plot_win(dod=dod)


def plot_window_sizes(wmin=10, wmax=25, step=5):
    for ii in range(wmin, wmax+1, step):
        print ii
        ObsData.winr = ii
        od = ObsData('flux', r'$S^{\rm Total}_{1.1} \ \ [{\rm Jy}]$')
        op = ObsPlot(od)
        op.plot_win()


def plot_all_obs():
    data = [
        # Flux density
        ObsData('flux', r'$S^{\rm Total}_{1.1} \ \ [{\rm Jy}]$'),
        ObsData('flux_40', r'$S_{1.1}(40^{\prime\prime}) \ \ [{\rm Jy}]$'),
        ObsData('flux_80', r'$S_{1.1}(80^{\prime\prime}) \ \ [{\rm Jy}]$'),
        ObsData('flux_120', r'$S_{1.1}(120^{\prime\prime}) \ \ [{\rm Jy}]$'),
        ObsData('fwhm_flux', r'$S^{\rm FWHM}_{1.1} \ \ [{\rm Jy}]$'),
        # Angular size
        ObsData('eqangled', r'$\theta^{\rm Total}_{\rm eq} \ \ [{\rm arcsec}]$'),
        ObsData('sangled', r'$\Omega^{\rm Total} \ \ [{\rm sq. \ arcsec}]$'),
        ObsData('fwhm_eqangled', r'$\theta^{\rm FWHM}_{\rm eq} \ \ [{\rm arcsec}]$'),
        ObsData('fwhm_sangled', r'$\Omega^{\rm FWHM} \ \ [{\rm sq. \ arcsec}]$'),
        ObsData('fwhm_sangled_ratio', r'$\Omega^{\rm FWHM} / \Omega^{\rm Total}$'),
        # NH3
        ObsData('nh3_tkin', r'$T_{\rm K} \ \ [{\rm K}]$'),
        ObsData('nh3_gbt_pk11', r'$T_{\rm pk}({\rm NH_3 \ (1,1)}) \ \ [{\rm K}]$'),
        ObsData('nh3_gbt_tau11', r'$\tau({\rm NH_3 \ (1,1)})$'),
        ObsData('nh3_gbt_fwhm', r'$\Delta v({\rm NH_3 \ (1,1)}) \ \ [{\rm km \ s^{-1}}]$'),
        # HCO+
        ObsData('mol_hco_tpk', r'$T_{\rm pk}({\rm HCO^+}) \ \ [{\rm K}]$'),
        ObsData('mol_hco_int', r'$I({\rm HCO^+}) \ \ [{\rm K \ km \ s^{-1}}]$'),
        ObsData('mol_hco_fwhm', r'$\Delta v({\rm HCO^+}) \ \ [{\rm km \ s^{-1}}]$'),
        ObsData('mol_hco_fwzi', r'${\rm HCO^+ \ FWZI} \ \ [{\rm km \ s^{-1}}]$'),
        # N2H+
        ObsData('mol_nnh_tpk', r'$T_{\rm pk}({\rm N_2H^+}) \ \ [{\rm K}]$'),
        ObsData('mol_nnh_int', r'$I({\rm N_2H^+}) \ \ [{\rm K \ km \ s^{-1}}]$'),
        ObsData('mol_nnh_fwhm', r'$\Delta v({\rm N_2H^+}) \ \ [{\rm km \ s^{-1}}]$'),
        ObsData('mol_nnh_fwzi', r'${\rm N_2H^+ \ FWZI} \ \ [{\rm km \ s^{-1}}]$'),
    ]
    for d in data:
        print ':: Plotting {0}'.format(d.col)
        op = ObsPlot(d)
        op.plot_frac()
        op.plot_tot()
        op.plot_win()


def plot_dist():
    od = DistData()
    op = DistPlot(od)
    op.plot_frac()


def plot_mass(use_fwhm=False):
    MassData.use_fwhm = use_fwhm
    od = MassData()
    op = MassPlot(od)
    op.plot_frac()


def plot_nfrac_mass():
    op = MassFracPlot(MassFracData())
    op.plot_frac()


def plot_lifetime(use_fwhm=False):
    LifeTimeData.use_fwhm = use_fwhm
    od = LifeTimeData()
    op = LifeTimePlot(od)
    op.plot_frac()


#####
# Old code to calculate fractions based on molecular detection.
#####

class SfFrac(object):
    df = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    nbins = 20.

    def __init__(self, col, df=None, cumsum=False, vmin=None, vmax=None):
        if df is not None:
            self.df = df
        assert col in self.df.columns
        assert self.df[col].shape[0] > 0
        self.col = col
        if vmin is not None:
            self.vmin = vmin
        else:
            self.vmin = self.df[col].min()
        if vmax is not None:
            self.vmax = vmax
        else:
            self.vmax = self.df[col].max()
        self.bins = np.linspace(self.vmin, self.vmax, self.nbins)
        self.step = (self.bins[1] - self.bins[0]) / 2.
        self.cumsum = cumsum

    def _get_hist(self):
        dets, edges = np.histogram(self.df[self.col], bins=self.bins)
        sf_dets, edges = np.histogram(self.df[self.df.sf_f == 1][self.col],
                                      bins=self.bins)
        if self.cumsum:
            self.dets = np.cumsum(dets)
            self.sf_dets = np.cumsum(sf_dets)
            self.frac = np.divide(self.sf_dets, 1.0 * self.dets)
        else:
            self.dets = dets
            self.sf_dets = sf_dets
            self.frac = np.divide(self.sf_dets, 1.0 * self.dets)

    def _make_fig(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

    def plot(self):
        self._get_hist()
        self._make_fig()
        self.ax.step(self.bins[:-1], self.frac, 'k-')
        #uncert = np.sqrt(self.sf_dets) / self.dets
        #self.ax.errorbar(self.bins[:-1] - self.step, self.frac, yerr=uncert,
        #        ecolor=0.4, fmt='.')
        self.ax.set_ylabel(r'$R_{\rm sf}$')


def plot_hco():
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    hco_dets = evo[evo.mol_hco_f.isin([1,3])]
    SfFrac.nbins = 30.
    # Cumulative
    sff1 = SfFrac('mol_hco_tpk', df=hco_dets, cumsum=True, vmin=0.1, vmax=3.0)
    sff1.plot()
    sff1.ax.set_xlabel(r'$T_{\rm pk}({\rm HCO^+})$')
    # Fraction per bin
    sff2 = SfFrac('mol_hco_tpk', df=hco_dets, cumsum=False, vmin=0.1, vmax=3.0)
    sff2.plot()
    sff2.ax.set_xlabel(r'$T_{\rm pk}({\rm HCO^+})$')
    return sff1, sff2


def plot_nnh():
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    nnh_dets = evo[evo.mol_nnh_f.isin([1,3])]
    SfFrac.nbins = 30.
    # Cumulative
    sff1 = SfFrac('mol_nnh_tpk', df=nnh_dets, cumsum=True, vmin=0.1, vmax=3.0)
    sff1.plot()
    sff1.ax.set_xlabel(r'$T_{\rm pk}({\rm N_2H^+})$')
    # Fraction per bin
    sff2 = SfFrac('mol_nnh_tpk', df=nnh_dets, cumsum=False, vmin=0.1, vmax=3.0)
    sff2.plot()
    sff2.ax.set_xlabel(r'$T_{\rm pk}({\rm N_2H^+})$')
    return sff1, sff2


def plot_flux():
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    SfFrac.nbins = 30.
    # Cumulative
    sff1 = SfFrac('flux_40', df=evo, cumsum=True, vmin=0.1, vmax=3.0)
    sff1.plot()
    sff1.ax.set_xlabel(r'$S_{\rm 1.1mm}(40^{\prime\prime})$')
    # Fraction per bin
    sff2 = SfFrac('flux_40', df=evo, cumsum=False, vmin=0.1, vmax=3.0)
    sff2.plot()
    sff2.ax.set_xlabel(r'$S_{\rm 1.1mm}(40^{\prime\prime})$')
    return sff1, sff2


