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
from besl import util


plt.rc('font', **{'size':14})


def get_data():
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    fwhm = catalog.read_cat('bgps_v210_fwhm').set_index('v210cnum')
    fwhm = fwhm.loc[:, 'npix':]
    df = evo.merge(fwhm, left_index=True, right_index=True)
    return df


class ObsData(object):
    winr = 25

    def __init__(self, col, xlabel, bool_proto=False):
        # labels
        self.col = col
        self.xlabel = xlabel
        self.sf_fcol = 'sf_frac_' + col
        self.sf_tcol = 'sf_tot_' + col
        self.sf_wcol = 'sf_win_' + col
        # use proto_prob or boolean probability
        if bool_proto:
            self.proto_col = 'is_proto'
        else:
            self.proto_col = 'proto_prob'
        # data
        self.df = self.get_df()
        self.mean_frac = self.df.is_proto.mean()
        self.cumix = self.df.index
        self.bins = self.df[col].values
        self.nbins = self.bins.shape[0]
        self.xmin, self.xmax = self.bins[0], self.bins[-1]
        # compute
        print ':: Computing {0}'.format(col)
        self.fvals, self.tvals, self.wvals = self.cumfrac(self.cumix)

    def get_df(self):
        df = get_data()
        df = df.query('10 < glon_peak < 65')
        df = df[df[self.col].notnull() & (df[self.col] > 0)].sort(self.col)
        stages, labels = dpdf_calc.evo_stages(bgps=df)
        df['is_proto'] = False
        df.loc[stages[1].index, 'is_proto'] = True
        for col in (self.sf_fcol, self.sf_tcol, self.sf_wcol):
            df[col] = np.nan
        return df

    def cumfrac(self, cumix):
        for nn, ii in enumerate(cumix):
            nproto = self.df.loc[:ii, self.proto_col].sum()
            self.df.loc[ii, self.sf_fcol] = nproto / nn
            self.df.loc[ii, self.sf_tcol] = nproto / self.nbins
            self.df.loc[ii, self.sf_wcol] = self.window(nn, cumix)
        return (self.df[self.sf_fcol].values,
                self.df[self.sf_tcol].values,
                self.df[self.sf_wcol].values)

    def window(self, nn, cumix):
        # repeats last value off edge
        if nn < self.winr:
            win = cumix.values[:self.winr+1]
        elif nn > self.nbins - self.winr:
            win = cumix.values[-self.winr:]
        else:
            win = cumix.values[nn-self.winr:nn+self.winr+1]
        return self.df.loc[win, self.proto_col].mean()


class ObsPlot(object):
    def __init__(self, od):
        self.od = od

    def make_fig(self):
        fig, ax = plt.subplots(figsize=(6, 2.5))
        ax.set_xscale('log')
        ax.set_xlabel(self.od.xlabel)
        ax.set_xlim(self.od.xmin, self.od.xmax)
        ax.set_ylim(0, 1.1)
        plt.subplots_adjust(bottom=0.22)
        return fig, ax

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

    def plot_win(self):
        fig, ax = self.make_fig()
        ax.set_ylabel(r'$R_{\rm proto}$')
        ax.annotate(r'${\rm Window} = ' + str(2 * self.od.winr) + '$',
                    xy=(0.725, 0.1), xycoords='axes fraction', fontsize=13)
        ax.hlines(1, self.od.xmin, self.od.xmax, linestyles='dashed',
                  colors='0.5')
        ax.hlines(self.od.mean_frac, self.od.xmin, self.od.xmax,
                  linestyles='dotted', colors='0.5')
        ax.vlines(self.od.bins, 1.02, 1.07, colors='black', linestyles='solid',
                  linewidth=0.1)
        ax.plot(self.od.bins, self.od.wvals, 'k-', drawstyle='steps')
        ax.plot(self.od.bins[:self.od.winr], self.od.wvals[:self.od.winr],
                color='0.5', linestyle='solid', drawstyle='steps')
        ax.plot(self.od.bins[-self.od.winr:], self.od.wvals[-self.od.winr:],
                color='0.5', linestyle='solid', drawstyle='steps')
        util.savefig('{0}_{1}'.format(self.od.sf_wcol, self.od.winr),
                     close=True)


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


class MargData(object):
    pass


class MargPlot(object):
    pass


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


