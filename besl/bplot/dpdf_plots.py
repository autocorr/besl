"""
==========
DPDF Plots
==========

Plotting routines for DPDFs and sampling.

"""

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from besl import catalog, units
from besl.dpdf_calc import mc_sampler_1d, gen_stages


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

def stages_hist(label, xlabel, bgps=[]):
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
    # evo stages
    stages = gen_stages(bgps=bgps, stages_group=2, label=label)
    # calculate lims and bins
    # TODO
    xmin = _np.nanmin([df[label].min() for df in stages])
    xmax = _np.nanmax([df[label].max() for df in stages])
    lbins = _np.logspace(_np.log10(xmin), _np.log10(xmax), 40)
    # plot settings
    kwargs_gen = {'bins': lbins, 'color': 'LightGray', 'histtype':
        'stepfilled', 'linestyle': 'solid', 'linewidth': 1}
    kwargs_labels = [
        {'label': 'Starless'},
        {'label': r'${\rm H_2O \ \  N}$'},
        {'label': r'${\rm IR \ \ Y}$'},
        {'label': r'${\rm H_2O \ \ Y}$'},
        {'label': r'${\rm EGO \ \ Y}$'},
        {'label': r'${\rm H\sc{II} \ \ Y}$'}]
    stages_labels = [ r'Starless', r'${\rm H_2O \ \  N}$', r'${\rm IR \ \ Y}$',
        r'${\rm H_2O \ \ Y}$', r'${\rm EGO \ \ Y}$', r'${\rm H\sc{II} \ \ Y}$']
    # create plot
    fig, axes = _plt.subplots(nrows=len(stages), ncols=1, sharex=True)
    for i, ax in enumerate(axes):
        # draw plots
        kwargs_hist = dict(kwargs_gen.items() + kwargs_labels[i].items())
        hist_heights, hist_edges = _np.histogram(stages[i][label], bins=lbins)
        ymax = _np.nanmax([_np.max(hist_heights), 1])
        df = stages[i][label]
        if len(df) > 0:
            ax.hist(df[_np.isfinite(df)].values, **kwargs_hist)
            ax.plot(df.median() * _np.ones(2), [0, 1.2 * ymax], 'k--')
        # plot attributes
        ax.set_xlim([10**(_np.log10(xmin) - 0.2), 10**(_np.log10(xmax) + 0.2)])
        ax.set_ylim([0, 1.1 * ymax])
        ax.locator_params(axis='y', tight=True, nbins=5)
        ax.set_xscale('log')
        #ax.legend(loc=1, frameon=False, numpoints=None, prop={'size':12})
        ax.annotate(stages_labels[i], xy=(0.875, 0.75), xycoords='axes fraction',
            fontsize=10)
    axes[-1].set_xlabel(xlabel)
    # save
    _plt.subplots_adjust(hspace=0.05)
    _plt.savefig('stages_hist_{}.pdf'.format(label))
    print '-- stages_hist_{}.pdf written'.format(label)
    return [fig, axes]

def write_all_stages_plots(bgps=[]):
    columns = [
        'flux',
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
        'rind_surf_area']
    dfs = [
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
        bgps[bgps.all_dML.notnull()]]
    labels = [
        r'$S_{1.1} \ \ [{\rm Jy}]$',
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
        r'${\rm Surface \ Area} \ \ [{\rm pc}]$']
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
    starless = bgps[(bgps.h2o_f == 0) & (bgps.corn_n == 0) & (bgps.ir_f == 0)]
    h2o_no = bgps[bgps.h2o_f == 0]
    ir_yes = bgps[bgps.ir_f == 1]
    h2o_yes = bgps[bgps.h2o_f == 1]
    hii_yes = bgps[bgps.corn_n > 0]
    ego_yes = bgps[bgps.ego_n > 0]
    df_list = [starless, h2o_no, ir_yes, h2o_yes, hii_yes, ego_yes]
    df_names = ['Starless', 'H2O No', 'IR Yes', 'H2O Yes', 'HII Yes', 'EGO Yes']
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


