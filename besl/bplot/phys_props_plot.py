#!/usr/bin/env python
# encoding: utf-8


from collections import namedtuple
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.optimize import leastsq
from matplotlib import colors, cm
from matplotlib import patheffects as PathEffects
from besl.dpdf_calc import evo_stages


class ScatterPanel(object):
    """
    """
    # TODO density contours
    marker_size = 3
    alpha = 0.5
    marker = '.'
    set_log = True

    def __init__(self, ax, df, cols, ax_labels=None, stage_label=None):
        for col in cols:
            assert col in df.columns
        self.ax = ax
        self.df = df
        self.cols = cols
        self.ax_labels = ax_labels
        self.stage_label = stage_label
        self.label_ypos = 0.85
        self.p0 = [1, 1]

    @staticmethod
    def _err_fun(p, y, x):
        return (y - (p[0] + x * p[1])) * (1 / y * 0.1)

    @staticmethod
    def _fit_fun(p, x):
        return p[0] + (x * p[1])

    @staticmethod
    def _power_law(p, x):
        return 10**p[0] * x**p[1]

    def get_fit(self):
        ndf = self.df[self.df[self.cols[0]].notnull() &
                      (self.df[self.cols[0]] > 0) &
                      self.df[self.cols[1]].notnull() &
                      (self.df[self.cols[1]] > 0)]
        x = ndf[self.cols[0]]
        y = ndf[self.cols[1]]
        self.p0[0] = y.median() / x.median()
        x = np.log10(x)
        y = np.log10(y)
        self.fit, self.fit_pars = leastsq(self._err_fun,
                                          self.p0[:], args=(y, x),
                                          maxfev=3000)

    def _draw_label(self, text, position=None):
        if position is None:
            position = (0.1, self.label_ypos)
            self.label_ypos -= 0.12
        txt = self.ax.annotate(text, xy=position,
                xycoords='axes fraction', fontsize=10)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2,
                foreground='w')])

    def set_stage_label(self):
        if self.stage_label is not None:
            self._draw_label(self.stage_label)

    def set_spear_label(self):
        dff = self.df[self.df[self.cols[0]].notnull() & self.df[self.cols[1]].notnull()]
        rho_spear, p_value = spearmanr(dff[self.cols])
        text = r'$\rho_{\rm sp} = ' + '{0:1.2f}$'.format(rho_spear)
        self._draw_label(text)

    def set_count_label(self):
        count = self.df[self.df[self.cols[0]].notnull() &
                        self.df[self.cols[1]].notnull()].shape[0]
        text = r'$N = {0}$'.format(count)
        self._draw_label(text)

    def set_pindex_label(self):
        text = r'$\beta = ' + '{0:1.2f}$'.format(self.fit[1])
        self._draw_label(text)

    def set_labels(self):
        if self.ax_labels is not None:
            self.ax.set_xlabel(self.ax_labels[0])
            self.ax.set_ylabel(self.ax_labels[1])

    def scatter(self, color=None):
        if color is None:
            color = 'black'
        self.ax.plot(self.df[self.cols[0]], self.df[self.cols[1]],
                marker=self.marker, markerfacecolor=color,
                markeredgecolor='None', markersize=self.marker_size,
                linestyle='None', alpha=self.alpha)
        if self.set_log:
            self.ax.set_xscale('log')
            self.ax.set_yscale('log')

    def plot_fit(self):
        try:
            self.get_fit()
            xmin = self.df[self.cols[0]].min()
            xmax = self.df[self.cols[0]].max()
            x = np.linspace(xmin, xmax, 1e3)
            self.ax.plot(x, self._power_law(self.fit, x), 'r-')
            self.set_pindex_label()
        except Exception:
            pass


class ScatterGrid(object):
    """
    """
    set_stage_label = True
    set_count_label = True
    set_spear = True
    cNorm = colors.Normalize(vmin=0, vmax=1)
    scalarmap = cm.ScalarMappable(norm=cNorm, cmap=cm.cubehelix)
    fig_dims = (2, 2)
    ncols = 4
    low_color_lim = 0.8  # 1 is white

    def __init__(self, ldf, cols, labels, ax_labels=None, lims=None, bw=False):
        assert len(ldf) > 0
        self.ldf = ldf
        self.cols = cols
        self.labels = labels
        self.ax_labels = ax_labels
        self.bw = bw
        # Plot related vars
        self.lims = lims
        self.panels = []
        self.rgb = [self.scalarmap.to_rgba(i) for i in np.linspace(0,
            self.low_color_lim, len(ldf))]

    def _get_axes(self):
        self.nrows = (len(self.ldf) // self.ncols) + 1
        figsize = (self.ncols * self.fig_dims[0],
                   self.nrows * self.fig_dims[1] + 0.2)
        self.fig, self.axes = plt.subplots(figsize=figsize, nrows=self.nrows,
                ncols=self.ncols, sharex=True, sharey=True)

    def _get_limits(self):
        if self.lims is None:
            self.lims = []
            self.lims.append(np.nanmin([df[self.cols[0]].min() for df in
                self.ldf]))
            self.lims.append(np.nanmax([df[self.cols[0]].max() for df in
                self.ldf]))
            self.lims.append(np.nanmin([df[self.cols[1]].min() for df in
                self.ldf]))
            self.lims.append(np.nanmax([df[self.cols[1]].max() for df in
                self.ldf]))

    def scatter(self):
        self._get_axes()
        self._get_limits()
        self.axes[0][0].set_xlim(self.lims[0], self.lims[1])
        self.axes[0][0].set_ylim(self.lims[2], self.lims[3])
        items = zip(self.ldf, self.rgb, self.labels, self.axes.flat)
        for df, color, label, ax in items:
            if self.bw:
                color = 'black'
            label = label.values()[0]
            panel = ScatterPanel(ax, df, self.cols, self.ax_labels,
                    stage_label=label)
            panel.scatter(color=color)
            if self.set_stage_label:
                panel.set_stage_label()
            if self.set_count_label:
                panel.set_count_label()
            if self.set_spear:
                panel.set_spear_label()
            #panel.plot_fit()
            if ax.is_first_col() & ax.is_last_row():
                panel.set_labels()
            self.panels.append(panel)
        self.axes.flatten()[len(self.ldf) - len(self.axes.flat)].axis('off')

    def save(self, png_only=False):
        if png_only:
            txt = '{0}_{1}.png'.format(*self.cols)
            plt.savefig(txt)
        else:
            for ext in ('pdf', 'eps', 'png'):
                txt = '{0}_{1}.{ext}'.format(*self.cols, ext=ext)
                plt.savefig(txt)


def plot_all_correlations(df, png_only=True):
    ldf, labels = evo_stages(df)
    plt.switch_backend('Agg')
    for ref1, ref2 in ColRefs.comb:
        sg = ScatterGrid(ldf, cols=[ref1.col, ref2.col],
                ax_labels=[ref1.ax_label, ref2.ax_label], labels=labels,
                bw=True)
        sg.scatter()
        sg.save(png_only=png_only)


class ColRefs(object):
    Ref = namedtuple('Ref', 'col ax_label')
    refs = [Ref('flux',
                r'$S_{270}$'),
            Ref('flux_40',
                r'$S_{270}(40^{\prime\prime})$'),
            Ref('h2o_gbt_int',
                r'$I({\rm H_2O})$'),
            Ref('h2o_gbt_vsp',
                r'$\Delta v({\rm H_2O})$'),
            Ref('h2o_gbt_num_lines',
                r'$N_{\rm pk}({\rm H_2O})$'),
            Ref('nh3_gbt_w11',
                r'$\Delta v_{\rm FWHM}({\rm NH_3})$'),
            Ref('mol_hco_tpk',
                r'$T_{\rm pk}({\rm HCO^+})$'),
            Ref('mol_hco_int',
                r'$I({\rm HCO^+})$'),
            Ref('mol_hco_fwhm',
                r'$\Delta v_{\rm FWHM}({\rm HCO^+})$'),
            Ref('mol_nnh_tpk',
                r'$T_{\rm pk}({\rm N_2H^+})$'),
            Ref('mol_nnh_int',
                r'$I({\rm N_2H^+})$'),
            Ref('nh3_tkin',
                r'$T_K({\rm NH_3})$'),
            Ref('ell_angle',
                r'$\theta$'),
            Ref('ell_sangle',
                r'$\Omega$'),
            Ref('dust_mass_surf',
                r'$\Sigma_{\rm d}$'),
            Ref('dML',
                r'$d_{\rm ML}$'),
            Ref('dust_mass',
                r'$M_{\rm d}$'),
            Ref('ell_area',
                r'$A$'),
            Ref('ell_radius',
                r'$R_{\rm eq}$'),
            Ref('vir_param',
                r'$\alpha_{\rm vir}$'),
            ]
    comb = combinations(refs, 2)


