#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from besl.catalog import read_cat


evo = read_cat('bgps_v210_evo').set_index('v210cnum')
hco_dets = evo[evo.mol_hco_f.isin([1,3])]
nnh_dets = evo[evo.mol_nnh_f.isin([1,3])]


class SfFrac(object):
    df = read_cat('bgps_v210_evo').set_index('v210cnum')
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
    SfFrac.nbins = 30.
    # Cumulative
    sff1 = SfFrac('mol_nnh_tpk', df=hco_dets, cumsum=True, vmin=0.1, vmax=3.0)
    sff1.plot()
    sff1.ax.set_xlabel(r'$T_{\rm pk}({\rm N_2H^+})$')
    # Fraction per bin
    sff2 = SfFrac('mol_nnh_tpk', df=hco_dets, cumsum=False, vmin=0.1, vmax=3.0)
    sff2.plot()
    sff2.ax.set_xlabel(r'$T_{\rm pk}({\rm N_2H^+})$')
    return sff1, sff2


def plot_flux():
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


