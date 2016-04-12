"""
=========================
Star Formation Flag Table
=========================

Star formation flag table based off of Rathborne et al. (2015) table of
correlation coefficients with colors in the bottom left.

"""
from __future__ import division
import itertools
from collections import OrderedDict
import numpy as np
import matplotlib as mpl
import matplotlib.patheffects as patheffects
from matplotlib import pyplot as plt
from besl import catalog
from besl import util


plt.rc('font', **{'size':10})


class Tile(object):
    lenf = 2  # data points width of square
    eps = 0.1  # spacing factor
    wid = lenf - 2 * eps
    cmap = plt.cm.YlGnBu
    norm = plt.Normalize(vmin=0, vmax=1)

    def __init__(self, f1, f2, ii, jj):
        self.f1 = f1
        self.f2 = f2
        self.uni_tot = f1.union(f2)
        self.int_tot = f1.intersection(f2)
        self.f1_tot = f1.sum()
        self.f2_tot = f2.sum()
        self.f1_frac = self.int_tot / self.f1_tot
        self.f2_frac = self.int_tot / self.f2_tot
        self.sq_pos = (self.lenf * ii + self.eps, self.lenf * jj + self.eps)

    def _get_pos(self, xx, yy):
        return (self.sq_pos[0] + xx * self.lenf, self.sq_pos[1] + yy * self.lenf)

    def add_patch(self, ax, cl):
        ax.add_patch(plt.Rectangle(self.sq_pos, self.wid, self.wid, facecolor=cl))

    def add_anno(self, ax, lstr, xx, yy):
        txt = ax.annotate(lstr, xy=self._get_pos(xx, yy), xycoords='data')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                                                     foreground='w')])

    def draw(self, ax):
        self.add_patch(ax, self.cmap(self.f1_frac))
        self.add_anno(ax, '{0}'.format(self.int_tot), 0.13, 0.51)
        self.add_anno(ax, '{0:1.3f}'.format(self.f1_frac), 0.13, 0.21)

    def draw_sum(self, ax):
        #self.add_patch(ax, '1.0')
        self.add_patch(ax, self.cmap(1.0))
        self.add_anno(ax, '{0}'.format(self.f1_tot), 0.13, 0.51)

    def draw_null(self, ax):
        ax.add_patch(plt.Rectangle(self.sq_pos, self.wid, self.wid,
                                   facecolor='0.5', hatch='/'))


class Flag(object):
    df = catalog.read_cat('bgps_v210_evo').set_index('v210cnum').query('10 < glon_peak < 65')

    def __init__(self, label, qstr):
        self.label = label
        self.qstr = qstr
        self.index = self.df.query(qstr).index

    def union(self, f):
        return self.index.union(f.index).shape[0]

    def intersection(self, f):
        return self.index.intersection(f.index).shape[0]

    def sum(self):
        return self.index.shape[0]


def get_all_flags(group=1):
    flags = []
    if group == 1:
        args = [
            (r'HG\ 0', 'hg70_eye_f == 0'),
            (r'HG\ 3', 'hg70_eye_f == 3'),
            (r'HG\ 1\&4', 'hg70_eye_f in [1,4]'),
            (r'HG\ 2', 'hg70_eye_f == 2'),
            # mid-IR
            (r'R08\ YSO', 'robit_f > 0'),
            (r'RMS', 'red_msx_f > 0'),
            (r'EGO', 'ego_n > 0'),
            # H2O
            (r'H_2O', 'h2o_f > 0'),
            # CH3OH
            (r'CH_3OH', 'ch3oh_f > 0'),
            # UCHII
            (r'UCHII', 'corn_n > 0'),
        ]
    elif group == 2:
        args = [
            (r'HG', 'hg70_eye_f in [1,2,4]'),
            # mid-IR
            (r'mid-IR', 'robit_f > 0 | red_msx_f > 0 | ego_n > 0'),
            # H2O
            (r'H_2O', 'h2o_f > 0'),
            # CH3OH
            (r'CH3_OH', 'ch3oh_f > 0'),
            # UCHII
            (r'UCHII', 'corn_n > 0'),
        ]
    else:
        raise ValueError('Invalid group: {0}.'.format(group))
    for label, qstr in args:
        flags.append(Flag(label, qstr))
    return flags


class GridPlot(object):
    figsize = (7, 7)

    def __init__(self, flags):
        self.flags = flags
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        self.ax.set_aspect('equal')

    def set_axes(self):
        labels = []
        for fl in self.flags:
            label = [''] + [r'${\rm ' + fl.label + r'}$']
            labels.extend(label)
        tvals = np.arange(0, len(labels) + 1)
        plt.xticks(tvals, labels, rotation=45, ha='right')
        plt.yticks(tvals, labels, rotation=45, ha='right', va='top')

    def hide_spines(self):
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.xaxis.set_ticks_position('none')
        self.ax.yaxis.set_ticks_position('none')

    def add_colorbar(self):
        box = self.ax.get_position()
        cax = plt.axes([
            box.x0 + box.width * 0.0065,
            box.y0 + box.height * 1.005,
            box.width * 0.3885,
            box.height * 0.05,
        ])
        mpl.colorbar.ColorbarBase(cax, norm=Tile.norm, cmap=Tile.cmap,
                                  ticks=[0, 0.25, 0.5, 0.75, 1],
                                  orientation='horizontal',
                                  ticklocation='top')

    def show_tiles(self):
        nn = len(self.flags)
        for ii, jj in itertools.product(range(nn), range(nn)):
            if ii == jj:
                Tile(self.flags[ii], self.flags[ii], ii, ii).draw_sum(self.ax)
            elif (ii < 4) & (jj < 4):
                Tile(self.flags[ii], self.flags[jj], ii, jj).draw_null(self.ax)
            else:
                Tile(self.flags[jj], self.flags[ii], jj, ii).draw(self.ax)

    def save(self, outname='intersect'):
        util.savefig(outname, close=True)


class RevGridPlot(GridPlot):
    def set_axes(self):
        labels = []
        for fl in self.flags:
            label = [''] + [r'${\rm ' + fl.label + r'}$']
            labels.extend(label)
        tvals = np.arange(0, len(labels) + 1)
        self.ax.xaxis.tick_top()
        plt.xticks(tvals, labels, rotation=45, ha='left')
        plt.yticks(tvals, [''] + labels[::-1][:-1], rotation=45, ha='right', va='top')

    def add_colorbar(self):
        box = self.ax.get_position()
        cax = plt.axes([
            box.x0 + box.width * 1.02,
            box.y0 - 0.0085,  # bottom right
            box.width * 0.025,
            box.height * 0.3775,
        ])
        mpl.colorbar.ColorbarBase(cax, norm=Tile.norm, cmap=Tile.cmap,
                                  ticks=[0, 0.25, 0.5, 0.75, 1],
                                  orientation='vertical',
                                  ticklocation='right')

    def show_tiles(self):
        nn = len(self.flags)
        for ii, jj in itertools.product(range(nn), range(nn)):
            iip = nn - ii - 1
            jjp = nn - jj - 1
            if ii == jj:
                Tile(self.flags[ii], self.flags[ii], ii, iip).draw_sum(self.ax)
            elif (ii < 4) & (jj < 4):
                Tile(self.flags[ii], self.flags[jj], ii, jjp).draw_null(self.ax)
            else:
                Tile(self.flags[jj], self.flags[ii], jj, iip).draw(self.ax)
        plt.subplots_adjust(bottom=0.05)


def make_plot():
    flags = get_all_flags(group=1)
    gp = RevGridPlot(flags)
    gp.set_axes()
    gp.hide_spines()
    gp.add_colorbar()
    gp.show_tiles()
    gp.save('rev_intersect')


