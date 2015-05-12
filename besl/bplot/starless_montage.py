#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
import aplpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from matplotlib import patheffects
from mpl_toolkits.axes_grid1.inset_locator import (zoomed_inset_axes, mark_inset)
from besl import catalog
from besl import image
from besl import util


plt.rc('font', **{'size':10})


class BgpsImg(object):
    label = r'$1.1 \ {\rm mm}$'
    beam = 0.009167
    vmin = -0.1
    vmax = 0.4
    co_sys = 'galactic'

    def __init__(self, cnum, lon, lat):
        self.img = self._get_image(cnum, lon, lat)

    def _get_image(self, cnum, lon, lat):
        img = image.get_bgps_img(cnum, 'map20')
        wcs = WCS(img[0].header)
        # crop to match HiGAL image
        #co_min = wcs.all_world2pix(lon - self.width / 2, lat - self.width / 2, 0)
        #co_max = wcs.all_world2pix(lon + self.width / 2, lat + self.width / 2, 0)
        #img[0].data = img[0].data[co_min[1].round():co_max[1].round(),
        #                          co_min[0].round():co_max[0].round()]
        return img

    def draw(self, fig, subplot):
        gc = aplpy.FITSFigure(self.img, figure=fig, convention='calabretta',
                              subplot=subplot)
        gc.show_colorscale(cmap='gist_gray', vmin=self.vmin, vmax=self.vmax,
                           stretch='linear')
        return gc


class HgImg(object):
    img_filen = '/home/svobodb/research/Data/HiGal/img_files/source/source_{0}_blue_svoboda.fits'
    label = r'$70 \ \mu {\rm m}$'
    beam = 0.0019444

    def __init__(self, cnum):
        self.img = fits.open(self.img_filen.format(cnum))
        self.vmin = np.percentile(self.img[0].data, 0.5)
        self.vmax = np.percentile(self.img[0].data, 90)

    def draw(self, fig, subplot):
        gc = aplpy.FITSFigure(self.img, figure=fig, convention='calabretta',
                              subplot=subplot)
        if self.vmax is not None:
            gc.show_colorscale(cmap='gist_gray', stretch='log',
                               vmin=self.vmin, vmax=self.vmax)
        else:
            gc.show_colorscale(cmap='gist_gray', stretch='linear')
        return gc


class ImgPair(object):
    cat = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')

    def __init__(self, cnum):
        self.cnum = cnum
        self.name = self.cat.loc[cnum, 'name'][7:]
        self.lon, self.lat = self.cat.loc[cnum, ['glon_cen', 'glat_cen']]
        self.mm = BgpsImg(cnum, self.lon, self.lat)
        self.hg = HgImg(cnum)


class PanelPair(object):
    label_fs = 10
    panhalfwid = 0.04166666

    def __init__(self, ii, cnum, fig):
        self.fig = fig
        self.data = ImgPair(cnum)
        self.mm_subplot = (2, 6, ii + 1)
        self.hg_subplot = (2, 6, ii + 6 + 1)

    def _recenter(self, gc):
        gc.recenter(self.data.lon, self.data.lat, radius=self.panhalfwid)

    def _panel_label(self, gc, label, xy=(0.08, 0.80)):
        txt = plt.annotate(label, xy=xy, xycoords='axes fraction',
                           fontsize=self.label_fs, weight='bold', color='black')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                                                     foreground='w')])

    def _add_cnum_label(self, gc):
        txt = plt.annotate(str(int(self.data.cnum)), xy=(0.6, 0.80),
                           xycoords='axes fraction', fontsize=self.label_fs,
                           weight='bold', color='black')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                foreground='w')])

    def _add_name_label(self, gc):
        txt = plt.annotate(self.data.name, xy=(0.08, 0.05),
                           xycoords='axes fraction', fontsize=self.label_fs,
                           weight='bold', color='black')
        txt.set_path_effects([patheffects.withStroke(linewidth=2,
                foreground='w')])

    def _add_beam(self, gc):
        gc.add_beam(facecolor='white', linewidth=2,
                    major=self.data.mm.beam,
                    minor=self.data.mm.beam, angle=0)

    def _hide_labels(self, gc):
        gc.axis_labels.hide_x()
        gc.axis_labels.hide_y()
        gc.tick_labels.hide_x()
        gc.tick_labels.hide_y()
        gc.ticks.hide_x()
        gc.ticks.hide_y()

    def _add_scalebar(self, gc):
        arcsec = 0.00833333
        gc.add_scalebar(arcsec)
        gc.scalebar.set_label(r'$30^{\prime\prime}$')
        gc.scalebar.set_linewidth(2)
        gc.scalebar.set_font_size(self.label_fs)
        gc.scalebar.set_color('w')

    def mm_panel(self):
        gc = self.data.mm.draw(self.fig, self.mm_subplot)
        self._recenter(gc)
        self._add_cnum_label(gc)
        self._hide_labels(gc)
        self._add_beam(gc)
        self._panel_label(gc, self.data.mm.label)
        self.mm_ax = gc

    def hg_panel(self):
        gc = self.data.hg.draw(self.fig, self.hg_subplot)
        self._hide_labels(gc)
        self._panel_label(gc, self.data.hg.label)
        self.mm_ax = gc

    def show(self):
        self.mm_panel()
        self.hg_panel()


class Montage(object):
    figsize = (7, 2.5)
    adjust_kwargs = {
        'bottom': 0.1,
        'top': 0.9,
        'left': 0.1,
        'right': 0.1,
        'hspace': 0.05,
        'wspace': 0.05,
    }

    def __init__(self, cnums):
        self.fig = plt.figure(figsize=self.figsize)
        self.cnums = cnums

    def show(self):
        ncols = 6
        nrows = 2
        for ii, cnum in enumerate(self.cnums):
            PanelPair(ii, cnum, self.fig).show()

    def save(self, outname='montage'):
        plt.tight_layout()
        util.savefig(outname, dpi=600, close=True)


def plot_starless_montage():
    sm = Montage([
        #5189,
        4029,
        4732,
        4119,
        5241,
        4788,
        4219,
    ])
    sm.show()
    sm.save()


