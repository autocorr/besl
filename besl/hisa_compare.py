"""
============
HISA Compare
============

Create images and analyze HI Self Absorption spectra.

"""

import os
import glob
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from astropy import wcs
from astropy.io import fits
from besl import catalog
from besl import util


omni = catalog.read_dpdf()
evo = catalog.read_cat('bgps_v210_evo').set_index('v210cnum')
velo = catalog.read_cat('bgps_v210_vel').set_index('v210cnum')
kdar = catalog.read_cat('bgps_v210_kdars').set_index('v210cnum')
cat = evo.loc[kdar.query('kdar == "F"').index].query('10 < glon_peak < 65')


class VgpsLib(object):
    """
    Container object for the VGPS survey data.

    Parameters
    ----------
    ext : string, default `".Tb.fits"`
        Extension of FITS image file.

    Attributes
    ----------
    abs_path : string
        Absolute path of data directory
    """
    abs_path = '/home/bsvoboda/data0/local_data/VGPS/'
    file_re = r'MOS_???'

    def __init__(self, ext='.Tb.fits'):
        self.ext = ext
        file_names = glob.glob(self.abs_path + self.file_re + ext)
        self.file_names = {int(os.path.basename(filen)[4:7]): filen
                           for filen in file_names}
        self.fields = sorted(self.file_names.keys())

    def get_hdu(self, field):
        hdu = fits.open(self.file_names[field])
        hdu[0].header['CUNIT3'] = 'm/s'
        return hdu

    def get_spectrum(self, hdu, lon, lat):
        """
        Extract spectrum at a given position of Lon and Lat.

        Parameters
        ----------
        hdu : `astropy.io.fits.HDU`
        lon, lat : float
            Galactic coordinates

        Returns
        -------
        spec : np.array
        """
        hdu_wcs = wcs.WCS(header=hdu[0].header)
        # glon, glat, vlsr, stokes, origin
        co = hdu_wcs.wcs_world2pix(lon, lat, 1, 1, 1)
        lon_pix, lat_pix, _, _ = [float(co_ii) for co_ii in co]
        spec = hdu[0].data[0, :, lat_pix, lon_pix]
        return spec

    def get_vaxis(self, hdu):
        """
        Get velocity axis from HDU cube.
        """
        vchan = hdu[0].data.shape[1]  # number of channels
        hdu_wcs = wcs.WCS(header=hdu[0].header)
        # glon, glat, vlsr, stokes, origin
        # Min velocity
        co = hdu_wcs.wcs_pix2world(1, 1, 1, 1, 1)
        _, _, min_v, _ = [float(co_ii) for co_ii in co]
        # Max velocity
        co = hdu_wcs.wcs_pix2world(1, 1, vchan, 1, 1)
        _, _, max_v, _ = [float(co_ii) for co_ii in co]
        return np.linspace(min_v, max_v, vchan)

    def get_border(self, hdu):
        """
        Get the coordinates of the rectangular edges of the image cube.

        Parameters
        ----------
        hdu : `astropy.io.fits.HDU`

        Returns
        -------
        coord : tuple
            (latmin, latmax, lonmin, lonmax)
        """
        shape = hdu[0].data.shape  # number of channels
        _, _, latn, lonn = shape
        hdu_wcs = wcs.WCS(header=hdu[0].header)
        # Bottom left corner
        co = hdu_wcs.wcs_pix2world(1, 1, 1, 1, 1)
        lonmax, latmin, _, _ = [float(co_ii) for co_ii in co]
        # Top right corner
        co = hdu_wcs.wcs_pix2world(lonn, latn, 1, 1, 1)
        lonmin, latmax, _, _ = [float(co_ii) for co_ii in co]
        return latmin, latmax, lonmin, lonmax


class HisaPlotter(object):
    """
    Draw HISA plots for a given dataset.

    Parameters
    ----------
    Lib : , default `VgpsLib`
        Data container object for the HISA HI data.

    lib_kwargs : dict, default `None`
        Keyword arguments to pass on to Lib

    Attributes
    ----------
    plot_dir : string
        Directory to save plots to.

    Usage
    -----
        >>> hp = HisaPlotter()
        >>> hp.draw_hisa_plots()
    """
    plot_dir = './plots/'

    def __init__(self, Lib=VgpsLib, lib_kwargs=None, verbose=True):
        if lib_kwargs is None:
            self.lib = Lib()
        else:
            self.lib = Lib(**lib_kwargs)
        self.verbose = verbose
        self.fig = plt.figure(figsize=(10, 5))
        self.ax = self.fig.add_subplot(111)

    def draw_hisa_plots(self):
        for field in self.lib.fields:
            if self.verbose:
                print '-- field: {0}'.format(field)
            hdu = self.lib.get_hdu(field)
            vaxis = self.lib.get_vaxis(hdu)
            latmin, latmax, lonmin, lonmax = self.lib.get_border(hdu)
            # FIXME this will break if goes through Galactic Center.
            field_cat = cat.query('@latmin < glat_peak < @latmax & '
                                      '@lonmin < glon_peak < @lonmax')
            for cnum in field_cat.index:
                lon = cat.loc[cnum, 'glon_peak']
                lat = cat.loc[cnum, 'glat_peak']
                vlsr = velo.loc[cnum, 'all_vlsr']
                spec = self.lib.get_spectrum(hdu, lon, lat)
                self.draw(cnum, spec, vaxis, vlsr)

    def draw(self, cnum, spec, vaxis, vlsr, to_kms=True):
        if to_kms:
            kms_factor = 1e-3
        else:
            kms_factor = 1
        # clip edges for noise
        spec = spec[10:-10]
        vaxis = vaxis[10:-10]
        self.ax.clear()
        self.ax.plot(vaxis * kms_factor, spec, 'k-')
        self.ax.vlines(vlsr, ymin=spec.min(), ymax=spec.max(),
                       linestyles='dashed')
        self.ax.vlines([vlsr - .1, vlsr + .1], ymin=spec.min(),
                       ymax=spec.max(), linestyles='dotted')
        self.ax.set_xlabel(r'$v_{\rm LSR} \ \ [{\rm km \ s^{-1}}]$')
        self.ax.set_ylabel(r'$T_{\rm b} \ \ [{\rm K}]$')
        self.ax.set_title(str(cnum))
        self.ax.grid()
        self.savefig(cnum)

    def savefig(self, cnum):
        if self.verbose:
            print '.',
        fbase = 'hisa_{0}'.format(cnum)
        for ext in ['.png', '.pdf', '.eps']:
            plt.savefig(self.plot_dir + fbase + ext, dpi=200)


