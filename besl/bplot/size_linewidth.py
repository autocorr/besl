"""
====================
Size Linewidth Plots
====================

Plotting routines for size linewidth relationships in the BGPS.

"""

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
import matplotlib.patheffects as PathEffects
from scipy.stats import spearmanr
from besl import catalog, dpdf_calc, util


def general_size_linewidth(stages, xcol, ycol, stages_labels, shape, colors,
    ax_labels):
    """
    """
    # TODO add doc
    # Plot settings
    xmin = _np.nanmin([df[xcol].min() for df in stages])
    xmax = _np.nanmax([df[xcol].max() for df in stages])
    ymin = _np.nanmin([df[ycol].min() for df in stages])
    ymax = _np.nanmax([df[ycol].max() for df in stages])
    # Begin plot
    error_kwargs = {'elinewidth': 0.5, 'ecolor': 'black', 'capsize': 0, 'fmt':
        'D', 'ms': 2.5}
    fig, axes = _plt.subplots(figsize=(shape[0] * 2, shape[1] * 2),
        nrows=shape[0], ncols=shape[1], sharex=True, sharey=True)
    for i, ax in enumerate(axes.flatten()):
        # Error bar plot
        x = stages[i][xcol].values
        y = stages[i][ycol].values
        xerr = stages[i][xcol].values / 20.
        yerr = stages[i][ycol + '_err'].values
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, color=colors[i], **error_kwargs)
        # Plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([10**(_np.log10(ymin) - 0.2), 10**(_np.log10(ymax) + 0.2)])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(ax_labels[0])
        ax.set_ylabel(ax_labels[1])
        ax.annotate(stages_labels[i], xy=(0.875, 0.75), xycoords='axes fraction',
            fontsize=10)
    _plt.subplots_adjust(hspace=0.05)
    _plt.savefig('size_linewidth_{0}_{1}.pdf'.format(ycol, _np.prod(shape)))
    return fig, ax

def spear_size_linewidth_four(stages):
    """
    """
    # TODO add doc
    ax_labels = [r'$R \ \ [{\rm pc}]$',
                 r'$\Delta v_{\rm HCO^+} \ \ [{\rm km \ s^{-1}}]$']
    stages_labels = [r'${\rm Starless}$',
                     r'${\rm H_2O \ \  N}$',
                     r'${\rm IR \ \ Y}$',
                     r'${\rm H_2O \ \ Y}$']
    colors = ['green', 'SlateBlue', 'red', 'DodgerBlue']
    xcol = 'avg_diam'
    ycol = 'hco_fwhm'
    # Plot limits
    stages = [df[(df[xcol].notnull()) & (df[ycol].notnull())] for df in stages]
    xmin = _np.nanmin([df[xcol].min() for df in stages])
    xmax = _np.nanmax([df[xcol].max() for df in stages])
    ymin = _np.nanmin([df[ycol].min() for df in stages])
    ymax = _np.nanmax([df[ycol].max() for df in stages])
    spears = [spearmanr(df[xcol].values, df[ycol].values) for df in stages]
    spears = [str(i[0])[:4] for i in spears]
    spears = [r'$\rho_{\rm spear} = ' + s + r'$' for s in spears]
    # Plot settings
    error_kwargs = {'elinewidth': 0.5, 'ecolor': 'black', 'capsize': 0, 'fmt':
        'D', 'ms': 2.5}
    # Begin plot
    fig, axes = _plt.subplots(figsize=(12, 4), nrows=1, ncols=4, sharex=True,
        sharey=True)
    for i, ax in enumerate(axes.flatten()):
        # Error bar plot
        x = stages[i][xcol].values
        y = stages[i][ycol].values
        xerr = stages[i][xcol].values * 0.1
        yerr = stages[i][ycol + '_err'].values
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, color=colors[i], **error_kwargs)
        linex = _np.linspace(0.01, 30, 100)
        liney = linex**(0.50) * 2.0
        ax.plot(linex, liney, 'k--', alpha=0.5)
        # Plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([10**(_np.log10(ymin) - 0.2), 10**(_np.log10(ymax) + 0.2)])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(ax_labels[0])
        if i == 0:
            ax.set_ylabel(ax_labels[1])
        stage_txt = ax.annotate(stages_labels[i], xy=(0.70, 0.90), xycoords='axes fraction',
            fontsize=10)
        spear_txt = ax.annotate(spears[i], xy=(0.65, 0.85), xycoords='axes fraction',
            fontsize=10)
        stage_txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
        spear_txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
    _plt.subplots_adjust(top=0.9, bottom=0.15, left=0.1, right=0.9, hspace=0.05,
        wspace=0.05)
    _plt.savefig('size_linewidth_{0}_{1}.pdf'.format('hco', '4panel'))
    return fig, axes

def spear_marginal_four(stages):
    """
    Plot histograms of the spearman rank for sampled draws from the DPDFs.
    """
    # Plot settings
    ax_labels = [r'$\rho_{\rm spear}$',
                 r'$f$']
    stages_labels = [r'${\rm Starless}$',
                     r'${\rm H_2O \ \  N}$',
                     r'${\rm IR \ \ Y}$',
                     r'${\rm H_2O \ \ Y}$']
    colors = ['green', 'SlateBlue', 'red', 'DodgerBlue']
    hist_kwargs = {'histtype': 'stepfilled', 'edgecolor': 'black', 'bins':
        _np.linspace(0, 1, 100)}
    xcol = 'avg_diam'
    ycol = 'hco_fwhm'
    # Calculate ranks
    good_kdars = ['T', 'F', 'N']
    stages = [df[(df[xcol].notnull()) & (df[ycol].notnull()) &
        ((df['neighbor_KDAR'].isin(good_kdars)) |
        (df['dpdf_KDAR'].isin(good_kdars)))] for df in stages]
    spears = [[], [], [], []]
    for i, stage in enumerate(stages):
        print i
        widths = stage[ycol].values
        # Draw distances
        radii_samples = dpdf_calc.gen_stage_area_samples(stage, nsample=1e4,
            radius=True, flatten=False) / 1e6
        # Calculate spearman rank for each draw
        for radii in radii_samples.T:
            spearman_rank = spearmanr(widths, radii)[0]
            spears[i].append(spearman_rank)
    # Begin plot
    fig, axes = _plt.subplots(figsize=(12, 1.5), nrows=1, ncols=4, sharex=True,
        sharey=True)
    for i, ax in enumerate(axes.flatten()):
        ax.hist(spears[i], facecolor=colors[i], **hist_kwargs)
        med_spear = _np.median(spears[i])
        ax.plot(med_spear, 40, 'Dk', markersize=5)
        spear_label = r'$\langle\rho_{\rm spear}\rangle_{1/2} = ' \
            + str(med_spear)[:4] + r'$'
        # Plot attributes
        if i == 0:
            ax.set_ylabel(ax_labels[1])
        ax.set_xlabel(ax_labels[0])
        ax.set_xticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([])
        stage_txt = ax.annotate(stages_labels[i], xy=(0.70, 0.75), xycoords='axes fraction',
            fontsize=10)
        spear_txt = ax.annotate(spear_label, xy=(0.55, 0.625), xycoords='axes fraction',
            fontsize=10)
        stage_txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
        spear_txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
    _plt.subplots_adjust(top=0.9, bottom=0.25, left=0.1, right=0.9, hspace=0.05,
        wspace=0.05)
    _plt.savefig('size_linewidth_spearman_{0}_{1}.pdf'.format('hco', '4panel'))
    return fig, axes

def write_stages_plot(bgps):
    """
    """
    bgps['h2o_vsp_err'] = 0.1
    stages_labels = [
        [r'Starless',
         r'Protostellar'],
        [r'Starless',
         r'${\rm IR \ \ Y, \ H_2O \ \ N}$',
         r'${\rm IR \ \ Y, \ H_2O \ \ Y}$'],
        [r'Starless',
         r'${\rm H_2O \ \  N}$',
         r'${\rm IR \ \ Y}$',
         r'${\rm H_2O \ \ Y}$',
         r'${\rm EGO \ \ Y}$',
         r'${\rm H\sc{II} \ \ Y}$']]
    stages_colors = [['green', 'red'],
              ['green', 'red', 'DodgerBlue'],
              ['green', 'SlateBlue', 'red', 'DodgerBlue', 'orange', 'blue']]
    stages_ax_labels = [[r'$R \ \ [{\rm pc}]$', r'$\Delta v_{\rm HCO^+} \ \ [{\rm km \ s^{-1}}]$'],
              [r'$R \ \ [{\rm pc}]$', r'$\Delta v_{\rm N_2H^+} \ \ [{\rm km \ s^{-1}}]$'],
              [r'$R \ \ [{\rm pc}]$', r'$\Delta v_{\rm H_2O} \ \ [{\rm km \ s^{-1}}]$']]
    shapes = [(1, 2), (1, 3), (2, 3)]
    for i in range(3):
        shape = shapes[i]
        stage_label = stages_labels[i]
        colors = stages_colors[i]
        ax_labels = stages_ax_labels[i]
        kwargs = {'shape': shape, 'stages_labels': stage_label, 'colors': colors,
                  'ax_labels': ax_labels}
        # HCO+ size linewidth
        hco_stages = dpdf_calc.gen_stages(bgps=bgps[bgps.hco_f.isin([1,3])],
            stages_group=i)
        fig, ax = general_size_linewidth(stages=hco_stages, xcol='avg_diam',
                ycol='hco_fwhm', **kwargs)
        # N2H+ size linewidth
        nnh_stages = dpdf_calc.gen_stages(bgps=bgps[bgps.nnh_f.isin([1,3])],
            stages_group=i)
        fig, ax = general_size_linewidth(stages=nnh_stages, xcol='avg_diam',
                ycol='nnh_fwhm', **kwargs)
        # H2O spread linewidth
        h2o_stages = dpdf_calc.gen_stages(bgps=bgps[bgps.h2o_gbt_f == 1],
            stages_group=i)
        fig, ax = general_size_linewidth(stages=h2o_stages, xcol='avg_diam',
                ycol='h2o_vsp', **kwargs)
    return


