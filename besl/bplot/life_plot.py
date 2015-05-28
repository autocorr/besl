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
MASS_BINS = np.logspace(MASS_LMIN, MASS_LMAX, MASS_NBIN)
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

    def calc_up_frac(self, counts):
        return [
            count[MASS_IMAX:].sum() / (counts[self.meth_ix][MASS_IMAX:]).sum()
            for count in counts
        ]


class GrowPops(Pops):
    grow_frac = 2.0

    def calc_counts(self):
        hist = lambda x, g : np.histogram(
            g * self.masses[x - 1], bins=self.bins)[0] / self.nsamples
        counts = [hist(ss.index, 1.0) for ii, ss in enumerate(self.stages)]
        counts[0] = hist(self.stages[0].index, self.grow_frac)
        return counts


class Panel(object):
    lo_lf = 2.5e4  # van der Walt (2005) class II lifetime
    hi_lf = 4.0e4
    md_lf = (lo_lf + hi_lf) / 2.
    xlim = (1e2, 1e5)
    ylim = (1e3, 3e6)
    bins = MASS_BINS[:MASS_IMAX]

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
        self.ax.fill_between([1e2, 4e2], 1e3, 1e7, facecolor='0.9',
                             edgecolor='none')

    def add_label(self, label):
        txt = self.ax.annotate(label, xy=(0.4, 0.85), xycoords='axes fraction')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                                                     foreground='w')])

    def add_hi_lim(self, up_frac, facecolor='black'):
        tau_lim = up_frac * self.md_lf
        self.ax.plot([10**4.1, 10**4.6], [tau_lim, tau_lim], 'k-', zorder=8)
        self.ax.plot([10**4.35], [tau_lim], marker='o',
                     markerfacecolor=facecolor, markersize=4, zorder=9)

    def plot_ff(self):
        ff = lambda x: 1.6558e7 * x**(-0.5)  # in yr, x in Msun
        self.ax.plot(self.xlim, [ff(x) for x in self.xlim], 'w-', zorder=6)
        self.ax.plot(self.xlim, [ff(x) for x in self.xlim], 'k--', zorder=7)

    def plot_life(self, vals, facecolor='0.3'):
        vals = vals[:MASS_IMAX]
        self.ax.fill_between(self.bins, self.lo_lf * vals, self.hi_lf * vals,
                             facecolor=facecolor, edgecolor='none', zorder=2)
        ol = 2
        self.ax.plot(self.bins, self.lo_lf * vals, 'k-', drawstyle='steps',
                     linewidth=ol, zorder=1)
        self.ax.plot(self.bins, self.hi_lf * vals, 'k-', drawstyle='steps',
                     linewidth=ol, zorder=1)
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


def life_plot():
    pops = Pops()
    gpops = GrowPops()
    plot = LifePlot(pops, gpops)
    plot.plot()
    plot.save()


