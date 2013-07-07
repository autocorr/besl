"""
====================
Size Linewidth Plots
====================

Plotting routines for size linewidth relationships in the BGPS.

"""

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
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

def spear_size_linewidth_four(stages, shape, colors, ax_labels):
    """
    """
    # TODO add doc
    ax_labels = [r'$\Delta v_{\rm HCO^+} \ \ [{\rm km \ s^{-1}}]$',
                 r'$R \ \ [{\rm pc}]$']
    stages_labels = [r'Starless',
                     r'${\rm H_2O \ \  N}$',
                     r'${\rm IR \ \ Y}$',
                     r'${\rm H_2O \ \ Y}$']
    colors = ['green', 'SlateBlue', 'red', 'DodgerBlue']
    xcol = 'hco_fwhm'
    ycol = 'avg_diam'
    # Plot settings
    xmin = _np.nanmin([df[xcol].min() for df in stages])
    xmax = _np.nanmax([df[xcol].max() for df in stages])
    ymin = _np.nanmin([df[ycol].min() for df in stages])
    ymax = _np.nanmax([df[ycol].max() for df in stages])
    # Begin plot
    error_kwargs = {'elinewidth': 0.5, 'ecolor': 'black', 'capsize': 0, 'fmt':
        'D', 'ms': 2.5}
    fig, axes = _plt.subplots(figsize=(8, 6), nrows=1, ncols=4, sharex=True,
        sharey=True)
    for i, ax in enumerate(axes.flatten()):
        # Error bar plot
        x = stages[i][xcol].values
        y = stages[i][ycol].values
        xerr = stages[i][xcol + '_err'].values
        yerr = stages[i][ycol].values / 20.
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, color=colors[i], **error_kwargs)
        # Plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([10**(_np.log10(ymin) - 0.2), 10**(_np.log10(ymax) + 0.2)])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(ax_labels[0])
        ax.set_ylabel(ax_labels[1])
        ax.annotate(stages_labels[i], xy=(0.70, 0.85), xycoords='axes fraction',
            fontsize=10)
    _plt.subplots_adjust(hspace=0.05)
    _plt.savefig('size_linewidth_{0}_{1}.pdf'.format('hco', '4panel'))
    return fig, ax

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

