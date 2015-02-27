#!/usr/bin/env python
"""
======================
Compare 70/24 um Plots
======================

Plot a comparison between the Hi-GAL and MIPSGAL images for a clump. For now,
Specifically for clump 5253.
"""
from __future__ import division
import aplpy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patheffects as PathEffects
from astropy.io import fits
from astropy.wcs import WCS
from besl import catalog
from besl import image


fdir = '/home/svobodb/research/temp/hg70_mips/'


class Panel(object):
    cnum = 5253  # v2.1.0
    glon = 30.62426  # deg
    glat = 0.5472881  # deg
    ra = 281.3195  # deg
    dec = -1.803885  # deg
    width = 0.05  # deg

    def __init__(self, fig):
        self.fig = fig
        self.cat = None
        self.arc = None
        self.gc = None

    def _get_cats(self):
        print ':: Reading in MIPSGAL cats'
        self.cat = catalog.read_cat('mipsgal_catalog_lclip')
        self.arc = catalog.read_cat('mipsgal_archive_lclip')
        llo = self.glon - self.width
        lhi = self.glon + self.width
        blo = self.glat - self.width
        bhi = self.glat + self.width
        self.cat = self.cat.query('@llo < l < @lhi & @blo < b < @bhi')
        self.arc = self.arc.query('@llo < l < @lhi & @blo < b < @bhi')

    def get_image(self):
        print ':: Reading in FITS'
        # FIXME uses hard path
        self.img = fits.open(fdir + self.filen)
        self.gc = aplpy.FITSFigure(self.img, figure=self.fig,
                                   subplot=self.subplot)

    def set_label(self):
        txt = plt.annotate(self.label, xy=(0.05, 0.055),
                           xycoords='axes fraction', color='black', fontsize=10)
        txt.set_path_effects([PathEffects.withStroke(linewidth=1.5,
                                                     foreground='w')])

    def show(self):
        self.get_image()
        self.gc.show_grayscale(**self.color_kwargs)
        self.show_markers()
        self.recenter()
        self.set_label()
        self.adjust()


class BgpsPanel(Panel):
    label = r'${\rm BGPS \ 1.1 \ mm}$'
    subplot = [0.13, 0.1, 0.35, 0.7]
    color_kwargs = {
        'vmin': -0.2,
        'vmax':  0.3,
        'stretch': 'linear',
        'smooth': None,
    }

    def show_markers(self):
        pass

    def get_image(self):
        print ':: Reading in FITS'
        self.img = image.get_bgps_img(self.cnum, exten='map20')
        self.gc = aplpy.FITSFigure(self.img, figure=self.fig,
                                   subplot=self.subplot,
                                   convention='calabretta')

    def show_rind(self):
        pass

    def recenter(self):
        self.gc.recenter(self.glon, self.glat, radius=self.width)

    def adjust(self):
        self.gc.tick_labels.set_font(size='x-small')
        self.gc.axis_labels.set_font(size='small')
        self.gc.tick_labels.set_xformat('dd.dd')
        self.gc.tick_labels.set_yformat('d.dd')
        self.gc.ticks.set_color('black')


class HigalPanel(Panel):
    filen = 'destripe_l030_blue_wgls_rcal.fits'
    label = r'${\rm Hi-GAL \ 70 \ \mu m}$'
    subplot = [0.13, 0.1, 0.35, 0.7]
    color_kwargs = {
        'vmin': 450,
        'vmax': 1300,
        'stretch': 'arcsinh',
        'smooth': None,
    }

    def show_markers(self):
        pass

    def show_rind(self):
        eps = 9e-1  # fudge factor for contouring
        rind = image.get_bgps_img(self.cnum, exten='labelmask')
        rind[0].data[np.isnan(rind[0].data)] = 0
        rind[0].data[rind[0].data != self.cnum] = 0
        rind[0].data[rind[0].data == self.cnum] = 1
        rind[0].data = rind[0].data.astype(float)
        self.gc.show_contour(rind, levels=[eps], colors='black',
                             convention='calabretta')
        self.gc.show_contour(rind, levels=[eps], colors='white',
                             linestyles='dashed', convention='calabretta')
        rind.close()

    def recenter(self):
        self.gc.recenter(self.glon, self.glat, radius=self.width)

    def adjust(self):
        self.gc.tick_labels.hide()
        self.gc.axis_labels.hide()
        self.gc.ticks.set_color('black')


class MipsgalPanel(Panel):
    filen = 'MG0310p005_024.fits'
    label = r'${\rm MIPSGAL \ 24 \ \mu m}$'
    subplot = [0.5, 0.1, 0.35, 0.7]
    color_kwargs = {
        'vmin': 25,
        'vmax': 41,
        'stretch': 'log',
        'smooth': None,
    }

    def show_markers(self):
        # markers
        self._get_cats()
        dx = -0.5 * 8e-3
        dy = -1.0 * 8e-3
        cdx = np.zeros(self.cat.ra.values.shape[0]) + dx * 0.7
        cdy = np.zeros(self.cat.ra.values.shape[0]) + dy * 0.7
        adx = np.zeros(self.arc.ra.values.shape[0]) + dx * 0.7
        ady = np.zeros(self.arc.ra.values.shape[0]) + dy * 0.7
        arrow = {'width': 2.5, 'head_width': 5, 'head_length': 3}
        self.gc.show_arrows(self.cat.ra.values - dx, self.cat.dec.values - dy,
                            cdx, cdy, edgecolor='black', facecolor='black',
                            linewidths=0.3, **arrow)
        self.gc.show_arrows(self.arc.ra.values - dx, self.arc.dec.values - dy,
                            cdx, cdy, edgecolor='white', facecolor='white',
                            linewidths=0.3, **arrow)

    def show_rind(self):
        pass

    def recenter(self):
        self.gc.recenter(self.ra, self.dec, radius=self.width / 2.35)

    def adjust(self):
        self.gc.tick_labels.hide()
        self.gc.axis_labels.hide()
        self.gc.ticks.hide()


class Plot(object):
    def save(self, outfilen='comp70'):
        for ext in ['png', 'pdf', 'eps']:
            name = outfilen + '.' + ext
            print ':: Saving to {0}'.format(name)
            self.fig.savefig(name, bbox_inches='tight', dpi=900)
        plt.close()


class TwoPlot(Plot):
    def __init__(self):
        self.fig = plt.figure(figsize=(8, 4))
        self.hp = HigalPanel(self.fig)
        self.mp = MipsgalPanel(self.fig)
        self.hp.show()
        self.mp.show()


class ThreePlot(Plot):
    def __init__(self):
        self.fig = plt.figure(figsize=(8, 2.7))
        BgpsPanel.subplot = [0.05, 0.1, 0.27, 0.8]
        HigalPanel.subplot = [0.35, 0.1, 0.27, 0.8]
        MipsgalPanel.subplot = [0.65, 0.1, 0.27, 0.8]
        self.bp = BgpsPanel(self.fig)
        self.hp = HigalPanel(self.fig)
        self.mp = MipsgalPanel(self.fig)
        self.bp.show()
        self.hp.show()
        self.mp.show()


def make_twoplot():
    tp = TwoPlot()
    tp.save(outfilen='comp70_twoplot')


def make_threeplot():
    tp = ThreePlot()
    tp.save(outfilen='comp70_threeplot')


