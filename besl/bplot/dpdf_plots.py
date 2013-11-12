"""
==========
DPDF Plots
==========

Plotting routines for DPDFs and sampling.

"""

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from matplotlib import colors, cm
from matplotlib.ticker import MaxNLocator
from besl import catalog, units
from besl.dpdf_calc import mc_sampler_1d, evo_stages, \
                           gen_stage_mass_samples, gen_stage_area_samples


def plot_dpdf_sampling(n=200):
    """
    Plot a Monte Carlo sampling of a DPDF
    """
    dpdf = catalog.read_dpdf()
    x = _np.arange(1000) * 20. + 20.
    y = dpdf[5].data[0]
    lims = [x.min(), x.max(), 0, 1]
    samples = mc_sampler_1d(x, y, lims=lims, nsample=1e4)
    fig = _plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(samples, bins=_np.linspace(lims[0], lims[1], 200), linewidth=2,
        histtype='stepfilled', normed=True, alpha=0.5, color='black')
    ax.plot(x, y / 20, 'k-', linewidth=2)
    ax.set_xlabel(r'$D \ \ [{\rm pc}]$')
    ax.set_ylabel(r'${\rm Probability}$')
    ax.set_ylim([0, y.max() * 1.1 / 20.])
    _plt.savefig('dpdf_test_sampling.pdf', format='pdf')
    return ax

def stages_hist(label, xlabel, bgps=None):
    """
    Create a histogram with the evolutionary stages overplotted.

    Parameters
    ----------
    label : string
        column name in bgps to plot
    xlabel : string
        string to write on X axis
    bgps : pd.DataFrame

    Returns
    -------
    fig : matplotlib.Figure
    ax : matplotlib.Axes
    """
    _plt.rc('font', **{'size':14})
    # evo stages
    stages, anno_labels = evo_stages(bgps=bgps, stages_group=2, label=label)
    # colors
    #colors = ['green', 'SlateBlue', 'red', 'DodgerBlue', 'Orange', 'magenta']
    cNorm = colors.Normalize(vmin=0, vmax=1)
    scalarmap = cm.ScalarMappable(norm=cNorm, cmap=cm.cubehelix)
    rgb = [scalarmap.to_rgba(i) for i in _np.linspace(0, 1, 7)]
    # calculate lims and bins
    xmin = _np.nanmin([df[label].min() for df in stages])
    xmax = _np.nanmax([df[label].max() for df in stages])
    lbins = _np.logspace(_np.log10(xmin), _np.log10(xmax), 40)
    # plot settings
    kwargs_gen = {'bins': lbins, 'color': 'LightGray', 'histtype':
        'stepfilled', 'linestyle': 'solid', 'linewidth': 1}
    stages_labels = [i.values()[0] for i in anno_labels]
    # create plot
    fig, axes = _plt.subplots(nrows=len(stages), ncols=1, sharex=True)
    for i, ax in enumerate(axes):
        # draw plots
        kwargs_hist = dict(kwargs_gen.items() + anno_labels[i].items())
        hist_heights, hist_edges = _np.histogram(stages[i][label], bins=lbins)
        ymax = _np.nanmax([_np.max(hist_heights), 1])
        df = stages[i][label]
        if len(df) > 0:
            # histogram
            ax.hist(df[_np.isfinite(df)].values, facecolor=rgb[i],
                **kwargs_hist)
            # vertical line for median
            ax.plot(df.median() * _np.ones(2), [0, 1.2 * ymax], 'w-')
            ax.plot(df.median() * _np.ones(2), [0, 1.2 * ymax], 'k--')
            # horizontal line for confidence interval
            ax.plot([_np.percentile(df.values, 10),
                     _np.percentile(df.values, 90)], _np.ones(2) * ymax, 'k:')
        # plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([0, 1.1 * ymax])
        ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax.locator_params(axis='y', tight=True, nbins=4)
        ax.set_xscale('log')
        #ax.legend(loc=1, frameon=False, numpoints=None, prop={'size':12})
        ax.annotate(stages_labels[i], xy=(0.775, 0.65), xycoords='axes fraction',
            fontsize=13)
        ax.annotate(df.shape[0], xy=(0.05, 0.65), xycoords='axes fraction',
            fontsize=13)
    axes[-1].set_xlabel(xlabel)
    # save
    _plt.subplots_adjust(hspace=0.05)
    _plt.savefig('stages_hist_{}.pdf'.format(label))
    print '-- stages_hist_{}.pdf written'.format(label)
    return [fig, axes]

def marginal_stages_hist(bgps=[], label='dust_mass', realiz=50, nsample=1e2):
    """
    Create a histogram of the marginalized clump dust mass with the
    evolutionary stages overplotted.

    Parameters
    ----------
    bgps : pd.DataFrame
    label : string
        Column to marginalize, valid options include 'dust_mass' and
        'rind_surf_area' for now.
    realiz : number
        Realizations to draw
    nsample : number
        Number of samples to draw for each source for each realization

    Returns
    -------
    fig : matplotlib.Figure
    ax : matplotlib.Axes
    """
    # evo stages
    stages, anno_labels = evo_stages(bgps=bgps, stages_group=2)
    if label == 'dust_mass':
        xlabel = r'$M_{\rm dust} \ \ [{\rm M_{\odot}}]$'
    elif label == 'rind_surf_area':
        xlabel = r'${\rm Area} \ \ [{\rm pc}^2]$'
    else:
        raise ValueError('Invalid label {0}.'.format(label))
    # calculate lims and bins
    xmin = _np.nanmin([df[label].min() for df in stages]) / 1.5
    xmax = _np.nanmax([df[label].max() for df in stages]) * 1.5
    lbins = _np.logspace(_np.log10(xmin), _np.log10(xmax), 100)
    # plot settings
    kwargs_gen = {'bins': lbins, 'color': 'Black', 'histtype':
        'step', 'linestyle': 'solid', 'linewidth': 1, 'alpha': 0.2}
    stages_labels = [i.values()[0] for i in anno_labels]
    # create plot
    fig, axes = _plt.subplots(nrows=len(stages), ncols=1, sharex=True)
    for i, ax in enumerate(axes):
        print '--', i
        # draw plots
        kwargs_hist = dict(kwargs_gen.items() + anno_labels[i].items())
        ymax = 0
        medians = []
        for j in range(realiz):
            print j
            if label == 'dust_mass':
                stage_samples = gen_stage_mass_samples(stages[i],
                    nsample=nsample)
            elif label == 'rind_surf_area':
                stage_samples = gen_stage_area_samples(stages[i],
                    nsample=nsample) / 1e6
            ax.hist(stage_samples, **kwargs_hist)
            hist_heights, hist_edges = _np.histogram(stage_samples, bins=lbins)
            top_bin = _np.nanmax([_np.max(hist_heights), 1])
            if top_bin > ymax:
                ymax = top_bin
            medians.append(_np.median(stage_samples))
        ax.plot(_np.median(medians) * _np.ones(2), [0, 1.2 * ymax], 'w-')
        ax.plot(_np.median(medians) * _np.ones(2), [0, 1.2 * ymax], 'k--')
        # plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([0, 1.1 * ymax])
        ax.locator_params(axis='y', tight=True, nbins=5)
        ax.set_yticklabels(labels=[])
        ax.set_xscale('log')
        #ax.legend(loc=1, frameon=False, numpoints=None, prop={'size':12})
        ax.annotate(stages_labels[i], xy=(0.875, 0.75), xycoords='axes fraction',
            fontsize=10)
        c_stage = stages[i]
        ax.annotate(c_stage[c_stage.all_dML.notnull()].shape[0], xy=(0.05,
            0.75), xycoords='axes fraction', fontsize=10)
    axes[-1].set_xlabel(xlabel)
    # save
    _plt.subplots_adjust(hspace=0.05)
    _plt.savefig('marg_stages_hist_{}.pdf'.format(label))
    print '-- marg_stages_hist_{}.pdf written'.format(label)
    return [fig, axes]

def write_all_stages_plots(bgps=[]):
    columns = [
        'flux',
        'flux_40',
        'hco_int',
        'hco_fwhm',
        'nnh_int',
        'nnh_fwhm',
        'h2o_tpk',
        'h2o_int',
        'h2o_vsp',
        'h2o_num_lines',
        'voffset_h2o_hco',
        'nh3_tkin',
        'rind_area',
        'all_dML',
        'dust_mass',
        'avg_diam',
        'rind_surf_area',
        'mass_surf_dens']
    dfs = [
        bgps,
        bgps,
        bgps[bgps.hco_f.isin([1,3])],
        bgps[bgps.hco_f.isin([1,3])],
        bgps[bgps.nnh_f.isin([1,3])],
        bgps[bgps.nnh_f.isin([1,3])],
        bgps[bgps.h2o_gbt_f > 0],
        bgps[bgps.h2o_gbt_f > 0],
        bgps[bgps.h2o_gbt_f > 0],
        bgps[bgps.h2o_gbt_f > 0],
        bgps[(bgps.h2o_gbt_f > 0) & (bgps.hco_f.isin([1,3]))],
        bgps[(bgps.nh3_pk11 / bgps.nh3_noise11 > 4) & (bgps.nh3_tkin < 3e2)],
        bgps,
        bgps[bgps.all_dML.notnull()],
        bgps[bgps.all_dML.notnull()],
        bgps[bgps.all_dML.notnull()],
        bgps[bgps.all_dML.notnull()],
        bgps]
    labels = [
        r'$S_{1.1} \ \ [{\rm Jy}]$',
        r'$S_{1.1}(40^{\prime\prime}) \ \ [{\rm Jy}]$',
        r'${\rm I(HCO^+) \ \ [K \ km \ s^{-1}}]$',
        r'${\rm HCO^+ \ FWHM \ \ [km \ s^{-1}}]$',
        r'${\rm I(N_2H^+) \ \ [K \ km \ s^{-1}}]$',
        r'${\rm N_2H^+ \ FWHM \ \ [km \ s^{-1}}]$',
        r'$T_{\rm pk}({\rm H_2O}) \ \ [{\rm K \ km \ s^{-1}}]$',
        r'${\rm I(H_2O) \ \ [K \ km \ s^{-1}}]$',
        r'${\rm H_2O} \ v_{\rm spread} \ \ [{\rm K \ km \ s^{-1}}]$',
        r'$N{\rm H_2O}$',
        r'$|v_{\rm HCO^+} - v_{\rm H_2O}| \ \ [{\rm km \ s^{-1}}]$',
        r'$T_{\rm K} \ \ [{\rm K}]$',
        r'${\rm Area} \ \ [{\rm arcsec^2}]$',
        r'${\rm dML} \ \ [{\rm kpc}]$',
        r'$M_{\rm dust} \ \ [M_{\odot}]$',
        r'${\rm Diameter} \ \ [{\rm pc}]$',
        r'${\rm Surface \ Area} \ \ [{\rm pc}]$',
        r'$\Sigma_{\rm H_2} \ \ [{\rm g \ cm^{-2}}]$']
    for col, label, df in zip(columns, labels, dfs):
        stages_hist(label=col, xlabel=label, bgps=df)
    return

def print_properties(bgps, out_filen='bgps_props.txt'):
    """
    Print an ASCII text file containing mean and median quantities from the
    BGPS all-matched catalog.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS all-matched catalog
    out_filen : string, default 'bgps_props.txt'
        Name of outfile
    """
    out_file = open(out_filen, 'w')
    df_list, anno_labels = evo_stages(bgps=bgps)
    df_names = [i.value()[0] for i in anno_labels]
    for i, df in enumerate(df_list):
        df['nnh/hco'] = df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           ((df.nnh_f == 1) | (df.nnh_f == 3))]['nnh_int'] / \
                        df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           ((df.nnh_f == 1) | (df.nnh_f == 3))]['hco_int']
        df['hco/nh3'] = df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['hco_int'] / \
                        df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['nh3_pk11']
        df['nnh/nh3'] = df[((df.nnh_f == 1) | (df.nnh_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['hco_int'] / \
                        df[((df.nnh_f == 1) | (df.nnh_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['nh3_pk11']
        # group numbers
        out_file.write('-- {}:\n'.format(df_names[i]))
        out_file.write('Number in group: {}\n'.format(
            df.shape[0]))
        out_file.write('Number with DPDFs: {}\n'.format(
            df[df.KDAR.isin(['N','F','T'])].shape[0]))
        out_file.write('Number with Tkin: {}\n'.format(
            df[(df.amm_f > 0) | (df.nh3_mult_n > 0)].shape[0]))
        out_file.write('Number with DPDF & Tkin: {}\n'.format(
            df[((df.amm_f > 0) | (df.nh3_mult_n > 0)) &
            (df.KDAR.isin(['N','F','T']))].shape[0]))
        for col in ['flux', 'flux_40', 'flux_80', 'flux_120']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[col].median(), df[col].shape[0]))
        for col in ['hco_tpk', 'hco_int', 'hco_fwhm']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[(df.hco_f == 1) | (df.hco_f == 3)][col].median(),
                df[(df.hco_f == 1) | (df.hco_f == 3)][col].shape[0]))
        for col in ['nnh_tpk', 'nnh_int', 'nnh_fwhm']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[(df.nnh_f == 1) | (df.nnh_f == 3)][col].median(),
                df[(df.nnh_f == 1) | (df.nnh_f == 3)][col].shape[0]))
        for col in ['nnh/hco', 'hco/nh3', 'nnh/nh3']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[col].median(), df[col].shape[0]))
        for col in ['h2o_tpk', 'h2o_int', 'h2o_vsp', 'h2o_num_lines']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[df.h2o_gbt_f > 0][col].median(), df[col].shape[0]))
        for col in ['nh3_tkin', 'nh3_pk11']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[df.nh3_pk11 / df.nh3_noise11 > 4][col].median(),
                df[df.nh3_pk11 / df.nh3_noise11 > 4][col].shape[0]))
        for col in ['dust_mass', 'avg_diam', 'rind_surf_area']:
            out_file.write('{0}: {1} ({2})\n'.format(col,
                df[df.KDAR.isin(['N','F','T'])][col].median(),
                df[df.KDAR.isin(['N','F','T'])][col].shape[0]))
    out_file.close()
    return


