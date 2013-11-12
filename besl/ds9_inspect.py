#!/usr/bin/env python
# encoding: utf-8

import os
import ds9
from astropy import wcs
from astropy.io import fits
from .catalog import read_cat
from .image import get_bgps_img


class Ds9FrameError(Exception):
    pass


class Inspector(object):
    rind_contour_color = 'yellow'
    bgps = read_cat('bgps_v210').set_index('v210cnum')
    v = 210

    def __init__(self, cnum):
        self.cnum = cnum
        self.glonp = self.bgps.ix[cnum, 'glon_peak']
        self.glatp = self.bgps.ix[cnum, 'glat_peak']
        self.glonc = self.bgps.ix[cnum, 'glon_cen']
        self.glatc = self.bgps.ix[cnum, 'glat_cen']
        self.d = ds9.ds9()
        self._set_init_prefs()

    def _set_init_prefs(self):
        self.d.set('prefs nancolor red')
        self.d.set('cmap grey')
        self.d.set('pan to {0} {1} wcs galactic'.format(
                   self.glonp, self.glatp))
        self.d.set('regions system wcs')
        self.d.set('regions sky galactic')
        self.d.set('regions skyformat degrees')

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

    def show_peak_cross(self):
        self.d.set('regions',
                   'galactic; cross point {0:.5f} {1:.5f}'.format(
                    self.glonp, self.glatp))

    def show_center_cross(self):
        self.d.set('regions',
                   'galactic; cross point {0:.5f} {1:.5f}'.format(
                   self.glonc, self.glatc))

    def zoom(self, zlevel):
        self.d.set('zoom to {0}'.format(zlevel))


class HiGalInspector(Inspector):
    """
    Visually inspect HiGal cutouts in DS9.
    """
    # Directories and file names
    root_dir = '/mnt/eld_data/HiGal/'
    data_dir = os.path.join(root_dir, 'dat_files')
    data_file = os.path.join(data_dir,
        'source_{cnum}_blue_svoboda_photall_err.dat')
    img_dir = os.path.join(root_dir, 'img_files')
    img_file = os.path.join(img_dir, '{img_type}',
        '{img_str}_{cnum}_blue_svoboda.fits')
    # Properties
    img_types = {'allder2': 'allder2_source',
                 'der1x': 'der1x_source',
                 'der1y': 'der1y_source',
                 'der2x': 'der2x_source',
                 'der2x45': 'der2x45_source',
                 'der2y' : 'der2y_source',
                 'der2y45': 'der2y45_source',
                 'source': 'source',
                 'mask_1': 'mask_1._source'}
    zlevel = 8
    max_scale = 1e4

    def __init__(self, cnum, img_type='source'):
        super(HiGalInspector, self).__init__(cnum)
        self.cnum = cnum
        if img_type not in self.img_types.keys():
            raise ValueError('Invalid img_type: {0}.'.format(img_type))
        self.img_type = img_type
        self.filen = self._format_img_infile()

    def _format_img_infile(self):
        return self.img_file.format(cnum=self.cnum, img_type=self.img_type,
                                    img_str=self.img_types[self.img_type])

    def _get_img(self):
        self.img = fits.open(self.filen)
        self.img_max = self.img[0].data.max()
        self.img_max = self.img[0].data.min()
        self.d.set_pyfits(self.img[0])

    def set_scale(self):
        self.d.set('scale linear')
        if self.img_max > self.max_scale:
            self.d.set('scale mode zmax')
        else:
            self.d.set('scale mode zscale')

    def source_update(self, cnum, img_type='source'):
        self.cnum = cnum
        self.img_type = img_type
        self.filen = self._format_img_infile()
        # FIXME
        # delete all frames
        # call __init__ base class
        # call self.view
        pass

    def view(self):
        self._get_img()
        self.show_source_contour()
        self.zoom(self.zlevel)
        self.set_scale()


