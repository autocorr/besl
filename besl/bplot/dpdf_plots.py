"""
==========
DPDF Plots
==========

Plotting routines for DPDFs and sampling.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib.ticker import MaxNLocator
from besl import catalog
from besl import util
from besl.dpdf_calc import (mc_sampler_1d, evo_stages,
                            gen_stage_mass_samples, gen_stage_area_samples)


plt.rc('font', **{'size':10})


def plot_dpdf_sampling(n=200):
    """
    Plot a Monte Carlo sampling of a DPDF
    """
    dpdf = catalog.read_dpdf()
    x = np.arange(1000) * 20. + 20.
    y = dpdf[5].data[0]
    lims = [x.min(), x.max(), 0, 1]
    samples = mc_sampler_1d(x, y, lims=lims, nsample=1e4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(samples, bins=np.linspace(lims[0], lims[1], 200), linewidth=2,
        histtype='stepfilled', normed=True, alpha=0.5, color='black')
    ax.plot(x, y / 20, 'k-', linewidth=2)
    ax.set_xlabel(r'$D \ \ [{\rm pc}]$')
    ax.set_ylabel(r'${\rm Probability}$')
    ax.set_ylim([0, y.max() * 1.1 / 20.])
    plt.savefig('dpdf_test_sampling.pdf', format='pdf')
    return ax


def plot_dpdf_sum(outname='evo_dpdfs'):
    """
    Create a histogram with the sum of the DPDFs per evolutionary stage.

    Returns
    -------
    fig : matplotlib.Figure
    ax : matplotlib.Axes
    """
    xdist = np.arange(1000) * 20. + 20.
    # read data
    pposts = catalog.read_pickle('ppv_dpdf_posteriors')
    dposts = catalog.read_pickle('dpdf_posteriors')
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    pevo = evo.loc[pposts.keys()]
    devo = evo.loc[dposts.keys()]
    pstages, anno_labels = evo_stages(bgps=pevo)
    dstages, anno_labels = evo_stages(bgps=devo)
    # plot
    fig, axes = plt.subplots(nrows=len(pstages), ncols=1, sharex=True,
                             figsize=(4, 3.5))
    for ii, ax in enumerate(axes.flat):
        if ii < 2:
            [spine.set_linewidth(2.0) for spine in ax.spines.itervalues()]
    for pdf, ddf, ax, alabel in zip(pstages, dstages, axes, anno_labels):
        # take posterior sum and normalize
        ax.tick_params(axis='y', which='major', labelsize=8)
        pydist = np.sum([pposts[ii] for ii in pdf.index], axis=0)
        pydist /= pydist.sum() * 20
        dydist = np.sum([dposts[ii] for ii in ddf.index], axis=0)
        dydist /= dydist.sum() * 20
        print pydist.sum() * 20, pydist.max()
        ax.plot(xdist * 1e-3, pydist * 1e3, 'k-', linewidth=2)
        ax.plot(xdist * 1e-3, dydist * 1e3, linestyle='dashed', linewidth=2,
            color='0.75')
        ax.fill_between(xdist * 1e-3, pydist * 1e3, color='0.5')
        ax.set_ylim([0, 0.35])
        ax.set_yticks([0, 0.1, 0.2, 0.3])
        # vertical line for median
        pmed = np.cumsum(pydist) * 20
        pmed = xdist[len(pmed[pmed < 0.5])] / 1e3
        ax.plot(pmed * np.ones(2), [0, 0.35], 'w-')
        ax.plot(pmed * np.ones(2), [0, 0.35], 'k--')
        # labels
        ax.annotate(alabel['label'], xy=(0.73, 0.55), xycoords='axes fraction',
                    fontsize=9)
        ax.annotate('$N = ' + str(pdf.shape[0]) + '$', xy=(0.45, 0.55),
                    xycoords='axes fraction', fontsize=9)
        ax.annotate('$(' + str(ddf.shape[0]) + ')$', xy=(0.625, 0.55),
                    xycoords='axes fraction', fontsize=9, color='0.5')
    axes[3].set_ylabel(r'${\rm Relative \ Probability}$')
    ax.set_xlabel(r'${\rm Heliocentric \ Distance \ \ [kpc]}$')
    plt.subplots_adjust(hspace=0.1, right=0.95, left=0.15, top=0.95,
                        bottom=0.125)
    util.savefig(outname)
    print '-- {0}.pdf written'.format(outname)
    return fig, axes


class StagesHist(object):
    totl_color = '0.50'
    dpdf_color = '0.25'
    histtype='stepfilled'

    def __init__(self, col, xlabel, bgps=None):
        """
        Create a histogram with the evolutionary stages overplotted.

        Parameters
        ----------
        label : string
            column name in bgps to plot
        xlabel : string
            string to write on X axis
        bgps : pd.DataFrame
        """
        if bgps is None:
            bgps = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
        else:
            assert bgps.index.name == 'v210cnum'
        self.col = col
        self.xlabel = xlabel
        self.bgps = bgps
        # other data
        self.post = catalog.read_pickle('ppv_dpdf_posteriors')
        self.dix = post.keys()
        self.stages, self.stage_labels = evo_stages(bgps=bgps, stages_group=2, label=label)
        self.nstages = len(self.stages)

    def create_figure(self):
        plt.rc('font', **{'size':14})
        self.fig, self.axes = plt.subplots(nrows=len(self.stages), ncols=1,
                                           sharex=True)

    def add_stage_labels(self):
        for ii, ax in enumerate(self.axes):
            ax.annotate(self.stages_labels[ii], xy=(0.775, 0.65),
                xycoords='axes fraction', fontsize=13)

    def set_stage_numbers(self):
        for ii, ax in enumerate(self.axes):
            num_str = r'$' + self.stages.shape[0] + r'$'
            ax.annotate(num_str, xy=(0.05, 0.65),
                xycoords='axes fraction', fontsize=13)

    def set_xlabel(self):
        ax = self.axes[-1]
        ax.set_xlabel(self.xlabel)

    def save(self, outname=None):
        """
        Parameters
        ----------
        outname : string, default None
            If left None, then it will be written as `stages_hist_{col}`
            where `col` is the column name taken from the DataFrame.
        """
        if outname is None:
            outname = 'stages_hist_{0}'.format(self.col)
        plt.savefig(outname + '.pdf')
        plt.savefig(outname + '.eps')
        plt.savefig(outname + '.png', dpi=300)
        print '-- {0} written'.format(outname)


class MargStagesHist(StagesHist):
    def __init__(self, fvals, dvals, col, xlabel):
        super(MargStagesHist, object).__init__(col, xlabel)
        self.fvals = fvals
        self.dvals = dvals

    def set_stage_numbers(self):
        for ii, ax in enumerate(self.axes):
            num_str = r'$' + self.stages.shape[0] + r'$'
            ax.annotate(num_str, xy=(0.05, 0.65),
                xycoords='axes fraction', fontsize=13)

    def create_figure(self):
        plt.rc('font', **{'size':14})
        self.fig, self.axes = plt.subplots(nrows=len(self.stages), ncols=1,
                                           sharex=True, sharey=True)

    def make_figure(self):
        self.create_figure()
        self.set_xlabel()
        self.save()


def stages_hist(label, xlabel, bgps=None, color=True, left_label=False):
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
    if bgps is None:
        bgps = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    else:
        assert bgps.index.name == 'v210cnum'
    # evo stages
    post = catalog.read_pickle('ppv_dpdf_posteriors')
    dix = post.keys()
    stages, anno_labels = evo_stages(bgps=bgps, stages_group=2, label=label)
    # colors
    if color:
        cNorm = colors.Normalize(vmin=0, vmax=1)
        scalarmap = cm.ScalarMappable(norm=cNorm, cmap=cm.cubehelix)
        rgb = [scalarmap.to_rgba(i) for i in np.linspace(0, 1, len(stages))]
    else:
        rgb = ['0.5' for _ in range(len(stages))]
    # calculate lims and bins
    xmin = np.nanmin([df[label].min() for df in stages])
    xmax = np.nanmax([df[label].max() for df in stages])
    lbins = np.logspace(np.log10(xmin), np.log10(xmax), 40)
    # plot settings
    kwargs_gen = {'bins': lbins, 'color': 'LightGray', 'histtype':
        'stepfilled', 'linestyle': 'solid', 'linewidth': 1}
    stages_labels = [i.values()[0] for i in anno_labels]
    # create plot
    fig, axes = plt.subplots(nrows=len(stages), ncols=1, sharex=True,
                             figsize=(4, 3.5))
    for i, ax in enumerate(axes):
        if i < 2:
            [spine.set_linewidth(2.0) for spine in ax.spines.itervalues()]
        # draw plots
        ax.tick_params(axis='y', which='major', labelsize=8)
        kwargs_hist = dict(kwargs_gen.items() + anno_labels[i].items())
        hist_heights, hist_edges = np.histogram(stages[i][label], bins=lbins)
        ymax = np.nanmax([np.max(hist_heights), 1])
        df = stages[i][label]
        if len(df) > 0:
            # histogram
            ax.hist(df[np.isfinite(df)].values, facecolor=rgb[i],
                **kwargs_hist)
            ax.hist(df[np.isfinite(df) & df.index.isin(dix)].values,
                bins=lbins, histtype='stepfilled', facecolor='0.25')
            # vertical line for median
            ax.plot(df.median() * np.ones(2), [0, 1.2 * ymax], 'w-')
            ax.plot(df.median() * np.ones(2), [0, 1.2 * ymax], 'k--')
            ax.plot(df[df.index.isin(dix)].median() * np.ones(2),
                [0, 1.2 * ymax], 'w-')
            ax.plot(df[df.index.isin(dix)].median() * np.ones(2),
                [0, 1.2 * ymax], 'k:')
            # horizontal line for confidence interval
            #ax.plot([np.percentile(df.values, 10),
            #         np.percentile(df.values, 90)], np.ones(2) * ymax, 'k:')
        # plot attributes
        ax.set_xlim([10**(np.log10(xmin) - 0.2), 10**(np.log10(xmax) + 0.2)])
        ax.set_ylim([0, 1.1 * ymax])
        ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax.locator_params(axis='y', tight=True, nbins=4)
        ax.set_xscale('log')
        if not left_label:
            slabel_xoffset = 0.73
        else:
            slabel_xoffset = 0.15
        ax.annotate(stages_labels[i], xy=(slabel_xoffset, 0.55),
                    xycoords='axes fraction', fontsize=9)
        ax.annotate(df.shape[0], xy=(0.05, 0.55), xycoords='axes fraction',
                    fontsize=9)
    axes[-1].set_xlabel(xlabel)
    axes[3].set_ylabel(r'$N$')
    # save
    plt.subplots_adjust(hspace=0.1, right=0.95, left=0.15, top=0.95,
                        bottom=0.125)
    util.savefig('stages_hist_{0}'.format(label))
    print '-- stages_hist_{}.pdf written'.format(label)
    return fig, axes


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
    xmin = np.nanmin([df[label].min() for df in stages]) / 1.5
    xmax = np.nanmax([df[label].max() for df in stages]) * 1.5
    lbins = np.logspace(np.log10(xmin), np.log10(xmax), 100)
    # plot settings
    kwargs_gen = {'bins': lbins, 'color': 'Black', 'histtype':
        'step', 'linestyle': 'solid', 'linewidth': 1, 'alpha': 0.2}
    stages_labels = [i.values()[0] for i in anno_labels]
    # create plot
    fig, axes = plt.subplots(nrows=len(stages), ncols=1, sharex=True)
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
            hist_heights, hist_edges = np.histogram(stage_samples, bins=lbins)
            top_bin = np.nanmax([np.max(hist_heights), 1])
            if top_bin > ymax:
                ymax = top_bin
            medians.append(np.median(stage_samples))
        ax.plot(np.median(medians) * np.ones(2), [0, 1.2 * ymax], 'w-')
        ax.plot(np.median(medians) * np.ones(2), [0, 1.2 * ymax], 'k--')
        # plot attributes
        ax.set_xlim([10**(np.log10(xmin) - 0.2), 10**(np.log10(xmax) + 0.2)])
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
    plt.subplots_adjust(hspace=0.05)
    plt.savefig('marg_stages_hist_{}.pdf'.format(label))
    print '-- marg_stages_hist_{}.pdf written'.format(label)
    return [fig, axes]


class PlotData(object):
    def __init__(self, label, xlabel, bgps=None, color=False, left_label=False):
        """
        label : str
            Dataframe column name
        xlabel : str
            LaTeX string to use as xlabel
        bgps : pd.DataFrame, Default None
        color : str, Default None
            Matplotlib color string to use
        left_label : bool, Default False
            Plot the label on the left instead of the right
        """
        self.label = label
        self.xlabel = xlabel
        if bgps is None:
            self.bgps = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
        else:
            assert bgps.index.name == 'v210cnum'
            self.bgps = bgps
        self.color = color
        self.left_label = left_label
        self.fig = None
        self.axes = None

    def plot(self):
        fig, axes = stages_hist(self.label, self.xlabel,
            bgps=self.bgps, color=self.color, left_label=self.left_label)
        self.fig = fig
        self.axes = axes

    def close(self):
        """ Close axes and figure instances. """
        plt.cla()
        plt.clf()
        plt.close()


def write_all_hists():
    """
    Read in evolutionary and fwhm data, then create histograms for all
    osberverable properties.
    """
    evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
    flux = evo.query('1e-2 < flux_40')
    fwhm = catalog.read_cat('bgps_v210_fwhm').set_index('v210cnum')
    fwhm = evo.merge(fwhm.loc[:, 'npix':], left_index=True, right_index=True)
    fwhm_res = fwhm[fwhm.fwhm_eqangled.notnull()]
    amm = evo.query('6 < nh3_tkin < 100')
    hco = evo.query('mol_hco_f in [1,2,3]')
    nnh = evo.query('mol_nnh_f in [1,2,3]')
    data = (
        # Flux data
        PlotData('flux', r'$S^{\rm Total}_{1.1} \ \ [{\rm Jy}]$', evo),
        PlotData('flux_40', r'$S_{1.1}(40^{\prime\prime}) \ \ [{\rm Jy}]$', flux),
        PlotData('flux_80', r'$S_{1.1}(80^{\prime\prime}) \ \ [{\rm Jy}]$', flux),
        PlotData('flux_120', r'$S_{1.1}(120^{\prime\prime}) \ \ [{\rm Jy}]$', flux),
        # FWHM properties
        PlotData('fwhm_flux', r'$S^{\rm FWHM}_{1.1} \ \ [{\rm Jy}]$', fwhm),
        PlotData('eqangled', r'$\theta^{\rm Total}_{\rm eq} \ \ [{\rm arcsec}]$', fwhm),
        PlotData('sangled', r'$\Omega^{\rm Total} \ \ [{\rm sq. \ arcsec}]$', fwhm),
        PlotData('fwhm_eqangled', r'$\theta^{\rm FWHM}_{\rm eq} \ \ [{\rm arcsec}]$', fwhm_res),
        PlotData('fwhm_sangled', r'$\Omega^{\rm FWHM} \ \ [{\rm sq. \ arcsec}]$', fwhm_res),
        PlotData('fwhm_sangled_ratio', r'$\Omega^{\rm FWHM} / \Omega^{\rm Total}$', fwhm_res, left_label=True),
        PlotData('fwhm_sangled_ratio_inv', r'$\Omega^{\rm Total} / \Omega^{\rm FWHM}$', fwhm_res),
        # NH3
        PlotData('nh3_tkin', r'$T_{\rm K} \ \ [{\rm K}]$', amm),
        PlotData('nh3_gbt_pk11', r'$T_{\rm pk}({\rm NH_3 \ (1,1)}) \ \ [{\rm K}]$', amm),
        PlotData('nh3_gbt_tau11', r'$\tau({\rm NH_3 \ (1,1)})$', amm),
        PlotData('nh3_gbt_fwhm', r'$\Delta v({\rm NH_3}) \ \ [{\rm km \ s^{-1}}]$', amm),
        # HCO+
        PlotData('mol_hco_tpk', r'$T_{\rm pk}({\rm HCO^+}) \ \ [{\rm K}]$', hco),
        PlotData('mol_hco_int', r'$I({\rm HCO^+}) \ \ [{\rm K \ km \ s^{-1}}]$', hco),
        PlotData('mol_hco_fwhm', r'$\Delta v({\rm HCO^+}) \ \ [{\rm km \ s^{-1}}]$', hco),
        PlotData('mol_hco_fwzi', r'${\rm HCO^+ \ FWZI} \ \ [{\rm km \ s^{-1}}]$', hco),
        # N2H+
        PlotData('mol_nnh_tpk', r'$T_{\rm pk}({\rm N_2H^+}) \ \ [{\rm K}]$', nnh),
        PlotData('mol_nnh_int', r'$I({\rm N_2H^+}) \ \ [{\rm K \ km \ s^{-1}}]$', nnh),
        PlotData('mol_nnh_fwhm', r'$\Delta v({\rm N_2H^+}) \ \ [{\rm km \ s^{-1}}]$', nnh),
        PlotData('mol_nnh_fwzi', r'${\rm N_2H^+ \ FWZI} \ \ [{\rm km \ s^{-1}}]$', nnh),
    )
    for d in data:
        d.plot()
        d.close()


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


