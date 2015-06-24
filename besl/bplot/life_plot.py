#!/usr/bin/env python
# encoding: utf-8
"""
======================
Lifetime Distributions
======================

Plot the lifetime of different categories by anchoring the relative ages to the
Class II methanol maser.

"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects
from besl import catalog
from besl import dpdf_calc
from besl import dpdf_mc
from besl import util
from besl.bplot import sf_frac


plt.rc('font', **{'size':10})
MASS_LMIN = 0.0 #2.3
MASS_LMAX = 5.0
MASS_NBIN = 400
MASS_BSIZ = (MASS_LMAX - MASS_LMIN) / MASS_NBIN
MASS_BINS = np.logspace(MASS_LMIN, MASS_LMAX, MASS_NBIN)
MASS_IMIN = int(MASS_NBIN * (3.0 - MASS_LMIN) / (MASS_LMAX - MASS_LMIN))
MASS_IMAX = int(MASS_NBIN * (4.0 - MASS_LMIN) / (MASS_LMAX - MASS_LMIN))


class Pops(object):
    nsamples = 1e4
    bins = MASS_BINS
    meth_ix = 5

    def __init__(self):
        self.cat = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
        self.stages, self.labels = dpdf_calc.evo_stages(self.cat)
        self.masses = dpdf_mc.MassSampler(self.cat, self.nsamples).draw()
        self.counts = self.calc_counts()
        self.counts.append(self.counts[0] + self.counts[1])
        self.frac_meth = self.counts[self.meth_ix] / self.counts[-1]
        self.up_frac = self.calc_up_frac(self.counts)
        self.rel_pop = self.calc_rel_pop(self.counts)
        self.rel_pop.append(self.rel_pop[0] + self.rel_pop[1])  # total
        self.labels.append({'label': r'${\rm All \ Clumps}$'})

    def calc_counts(self):
        return [
            np.histogram(self.masses[stage.index - 1], bins=self.bins)[0]
                / self.nsamples
            for stage in self.stages
        ]

    def calc_rel_pop(self, counts):
        return [count / counts[self.meth_ix] for count in counts]

    def calc_up_frac(self, counts, low_bin=MASS_IMAX):
        return [
            count[low_bin:].sum() / (counts[self.meth_ix][low_bin:]).sum()
            for count in counts
        ]

    def calc_up_mass(self, low_mass):
        low_bin = np.argmin(np.abs(self.bins - low_mass))
        return self.calc_up_frac(self.counts, low_bin=low_bin)

    def calc_at_mass(self, at_mass):
        at_bin = np.argmin(np.abs(self.bins - at_mass))
        return [pop[at_bin] for pop in self.rel_pop]

    def calc_up_count(self, low_mass):
        low_bin = np.argmin(np.abs(self.bins - low_mass))
        return [count[low_bin:].sum() for count in self.counts]


class GrowPops(Pops):
    grow_frac = 2.0

    def calc_counts(self):
        hist = lambda x, g : np.histogram(
            g * self.masses[x - 1], bins=self.bins)[0] / self.nsamples
        counts = [hist(ss.index, 1.0) for ii, ss in enumerate(self.stages)]
        counts[0] = hist(self.stages[0].index, self.grow_frac)
        return counts


class Panel(object):
    # lifetimes from van der Walt (2005) scaled by MW SFR
    lo_lf = 3.5e4 * (4 / 2.3)
    hi_lf = 3.5e4 * (4 / 1.5)
    md_lf = (lo_lf + hi_lf) / 2.
    xlim = (1e2, 1e5)
    ylim = (3e3, 1e7)
    bins = MASS_BINS[MASS_IMIN:MASS_IMAX]

    def __init__(self, ax, vals, label, up_frac):
        self.ax = ax
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.grid()
        self.show_incomplete()
        self.add_label(label)
        self.add_hi_lim(up_frac)
        self.plot_life(vals)
        self.plot_ff()

    def show_incomplete(self):
        self.ax.fill_between([1e2, 4e2], 1e3, 1e8, facecolor='0.9',
                             edgecolor='none')

    def add_label(self, label):
        txt = self.ax.annotate(label, xy=(0.4, 0.85), xycoords='axes fraction')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                                                     foreground='w')])

    def add_hi_lim(self, up_frac, facecolor='black'):
        tau_lim = up_frac * self.md_lf
        self.ax.plot([10**4.35, 10**4.35],
                     [up_frac * self.lo_lf, up_frac * self.hi_lf], 'k-',
                     zorder=8)
        self.ax.plot([10**4.1, 10**4.6], [tau_lim, tau_lim], 'k-', zorder=8)
        self.ax.plot([10**4.35], [tau_lim], marker='o',
                     markerfacecolor=facecolor, markersize=4, zorder=9)

    def plot_ff(self):
        ff = lambda x: 1.6558e7 * x**(-0.5)  # in yr, x in Msun
        self.ax.plot(self.xlim, [ff(x) for x in self.xlim], 'w-', zorder=6)
        self.ax.plot(self.xlim, [ff(x) for x in self.xlim], 'k--', zorder=7)

    def plot_life(self, vals, facecolor='0.3'):
        vals = vals[MASS_IMIN:MASS_IMAX]
        print ':: Mean slope: ', (np.mean(np.log10(vals[1:]) -
                                          np.log10(vals[:-1])) / MASS_BSIZ)
        self.ax.fill_between(self.bins, self.lo_lf * vals, self.hi_lf * vals,
                             facecolor=facecolor, edgecolor='none', zorder=2)
        ol = 2
        self.ax.plot(self.bins, self.lo_lf * vals, 'k-', drawstyle='steps',
                     linewidth=ol, zorder=1)
        self.ax.plot(self.bins, self.hi_lf * vals, 'k-', drawstyle='steps',
                     linewidth=ol, zorder=1)
        #self.ax.plot(self.bins, self.md_lf * vals, 'r-', drawstyle='steps',
        #             linewidth=1.2, zorder=3)
        # low value cap
        self.ax.plot([self.bins[0], self.bins[0]],
                     [self.lo_lf * vals[0], self.hi_lf * vals[0]], 'k-',
                     linewidth=ol, zorder=1)
        # high value cap
        self.ax.plot([self.bins[-1], self.bins[-1]],
                     [self.lo_lf * vals[-1], self.hi_lf * vals[-1]], 'k-',
                     linewidth=ol, zorder=1)


class LifePlot(object):
    def __init__(self, pops, gpops=None):
        self.pops = pops
        self.gpops = gpops
        self.fig, self.axes = plt.subplots(nrows=2, ncols=4, figsize=(8, 4),
                                           sharex=True, sharey=True)
        plt.subplots_adjust(left=0.07, right=0.96, top=0.95, bottom=0.13,
                            hspace=0.15, wspace=0.15)

    def add_axes_labels(self):
        ax = self.axes[1, 0]
        ax.set_xlabel(r'$M_{\rm H_2}^{\rm Total} \ \ [{\rm M_\odot}]$')
        ax.set_ylabel(r'$\tau_{\rm phase} \ \ [{\rm yr}]$')

    def plot(self):
        ix = [0, 1, 7, 0, 2, 3, 4, 6]
        for ii, ax in enumerate(self.axes.flat):
            jj = ix[ii]
            if ii <= 2:
                [spine.set_linewidth(2) for spine in ax.spines.itervalues()]
            if ii == 3:
                ax.axis('off')
                continue
            pn = Panel(ax, self.pops.rel_pop[jj],
                       self.pops.labels[jj]['label'],
                       self.pops.up_frac[jj])
            if (ii in [0, 2]) & (self.gpops is not None):
                color = '0.75'
                pn.add_hi_lim(self.gpops.up_frac[jj], facecolor=color)
                pn.plot_life(self.gpops.rel_pop[jj], facecolor=color)
        self.add_axes_labels()

    def save(self, outname='life_grid'):
        util.savefig(outname, close=True)


class TwoLifePlot(LifePlot):
    def __init__(self, pops, gpops=None):
        self.pops = pops
        self.gpops = gpops
        self.fig, self.axes = plt.subplots(nrows=1, ncols=2, figsize=(4, 2.2),
                                           sharex=True, sharey=True)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.22, top=0.95,
                            hspace=0.05, wspace=0.07)

    def add_axes_labels(self):
        ax = self.axes[0]
        ax.set_xlabel(r'$M_{\rm H_2}^{\rm Total} \ \ [{\rm M_\odot}]$')
        ax.set_ylabel(r'$\tau_{\rm phase} \ \ [{\rm yr}]$')

    def plot(self):
        ix = [0, 7]
        for ii, ax in enumerate(self.axes.flat):
            jj = ix[ii]
            if ii <= 2:
                [spine.set_linewidth(2) for spine in ax.spines.itervalues()]
            pn = Panel(ax, self.pops.rel_pop[jj],
                       self.pops.labels[jj]['label'],
                       self.pops.up_frac[jj])
            if (ii in [0, 1]) & (self.gpops is not None):
                color = '0.75'
                pn.add_hi_lim(self.gpops.up_frac[jj], facecolor=color)
                pn.plot_life(self.gpops.rel_pop[jj], facecolor=color)
        self.add_axes_labels()
        for ax in self.axes.flat:
            ax.set_xlim([4e2, 1e5])

    def save(self, outname='life_grid_two'):
        util.savefig(outname, close=True)


def life_plot():
    pops = Pops()
    gpops = GrowPops()
    plot = LifePlot(pops, gpops)
    plot.plot()
    plot.save()


def life_plot_two():
    pops = Pops()
    gpops = GrowPops()
    plot = TwoLifePlot(pops, gpops)
    plot.plot()
    plot.save()


def print_up_lifes():
    pops = Pops()
    gpops = GrowPops()
    masses = np.logspace(2, 4, 50)
    print_str = '-- M: {0:6.0f}, T_up: {1:6.0f} T_at: {2:6.0f}'
    print '\n\n:: No Growth'
    for mm in masses:
        tup = 7.4e4 * pops.calc_up_mass(mm)[0]
        tat = 7.4e4 * pops.calc_at_mass(mm)[0]
        print print_str.format(mm, tup, tat)
    print '\n\n:: With Growth'
    for mm in masses:
        tup = 7.4e4 * gpops.calc_up_mass(mm)[0]
        tat = 7.4e4 * gpops.calc_at_mass(mm)[0]
        print print_str.format(mm, tup, tat)
    print '\n\n:: Counts'
    for cc in pops.calc_up_count(650):
        print '-- C: {0:4.2f}'.format(cc)


def plot_growth_cdf():
    # draw masses
    cat = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    stages, labels = dpdf_calc.evo_stages(cat)
    masses = dpdf_mc.MassSampler(cat, 1e2).draw()
    # transform masses
    def bondi(m):
        prefac = 3 * np.pi / (8 * np.sqrt(2)) * 0.03**0.25
        t_clump = 7.5e5 * (m / 1e3)**(-0.5)
        t_cloud = 2e6
        tau = np.minimum(t_clump, t_cloud) / t_cloud
        return m * (1 - prefac * tau)**(-4.0)
    scc_obs = masses[stages[0].index - 1].flatten()
    pro_obs = masses[stages[1].index - 1].flatten()
    scc_bon = bondi(scc_obs)
    # plotting
    bins = np.logspace(1, 5, 400)
    fig, ax = plt.subplots(figsize=(4.0, 2.0))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.25)
    hist_args = {'bins': bins, 'cumulative': True, 'normed': True,
                 'histtype': 'step'}
    ax.plot([bins[0], bins[-1]], [1, 1], color='0.75', linestyle='dashed',
            zorder=0)
    ax.hist(scc_obs, color='cyan', zorder=2, **hist_args)
    ax.hist(2.0 * scc_obs, color='0.75', **hist_args)
    ax.hist(2.5 * scc_obs, color='0.50', **hist_args)
    ax.hist(3.0 * scc_obs, color='0.25', **hist_args)
    ax.hist(scc_bon, color='magenta', **hist_args)
    ax.hist(pro_obs, color='red', linestyle='dashed', **hist_args)
    ax.set_ylim([0, 1.1])
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_{\rm H_2}^{\rm Total} \ \ [{\rm M_\odot}]$')
    ax.set_ylabel(r'${\rm CDF}$')
    util.savefig('growth_cdf', close=True)


