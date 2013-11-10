#!/usr/bin/env python
# encoding: utf-8

import os
import ds9
from astropy import wcs
from astropy.io import fits
from .image import get_bgps_img


class Ds9FrameError(Exception):
    pass


class Inspector(object):
    rind_contour_color = 'yellow'
    v = 210

    def __init__(self, cnum):
        self.cnum = cnum
        self.d = ds9.ds9()

    def _get_frame_ids(self):
        frame_ids = self.d.get('frame all')
        frame_ids = [int(i) for i in frame_ids.split()]
        return frame_ids

    def _get_next_frame_id(self):
        frame_ids = self._get_frame_ids()
        next_frame_id = max(frame_ids) + 1
        if next_frame_id > 38:
            raise Ds9FrameError
        return next_frame_id

    def _get_rind(self):
        rind = get_bgps_img(self.cnum, exten='labelmask', v=self.v)
        rind[0].data[rind[0].data != self.cnum] = 0
        return rind

    def _get_flux(self):
        flux = get_bgps_img(self.cnum, exten='map20', v=self.v)
        return flux

    def show_source_contour(self):
        rind = self._get_rind()
        fid = self._get_next_frame_id()
        self.d.set('frame {0}'.format(fid))
        self.d.set_pyfits(rind)
        self.d.set('contour limits 0 1')
        self.d.set('contour nlevels 1')
        self.d.set('contour smooth 1')
        self.d.set('contour method block')
        self.d.set('contour')
        self.d.set('contour copy')
        self.d.set('frame 1')
        self.d.set('contour paste {0}'.format(self.rind_contour_color))
        self.d.set('frame delete {0}'.format(fid))

    def show_flux_contour(self, clevels):
        fid = self._get_next_frame_id()
        self.d.set('frame {0}'.format(fid))
        flux_img = self._get_flux()
        self.d.set_pyfits(flux_img)
        self.d.set('contour nvelels {0}'.format(len(clevels)))
        self.d.set('contour levels {' + ' '.join([str(i) for i in clevels]) + '}')
        self.d.set('contour smooth 2')
        self.d.set('contour')
        self.d.set('contour copy')
        self.d.set('frame 1')
        self.d.set('contour paste cyan')
        self.d.set('frame delete {0}'.format(fid))

    def zoom(self, xcen, ycen, xwidth, ywidth):
        pass


class HiGalInspector(Inspector):
    """
    Visually inspect HiGal cutouts in DS9.
    """
    # Directories and file names
    root = '/mnt/eld_data/HiGal/'
    data = os.path.join(self.root, 'dat_files')
    data_file = os.path.join(self.data,
        'source_{cnum}_blue_svoboda_photall_err.dat')
    img = os.path.join(self.root, 'img_files')
    img_file = os.path.join(self.img, '{img_type}',
        '{img_type}_source_{cnum}_blue_svoboda.fits')
    # Properties
    img_types = ['allder2', 'der1x', 'der1y', 'der2x', 'der2x45', 'der2y',
                 'der2y45', 'mask_1.']

    def __init__(self, cnum):
        super(HigalInspector, self).__init__(cnum)
        self.cnum = cnum


    def get_img(self, cnum, img_type):
        self.filen = self.img_file.format(cnum=cnum, img_type=img_type)
        self.img = fits.open(self.filen)
        self.d.set_pyfits(self.img)

    def view(self, cnum, img_type):
        self.show_rind(cnum)
        self.get_img(cnum, img_type)


