#!/usr/bin/env python
# encoding: utf-8

import os
import pandas as pd
import pyds9
from astropy.io import fits
from .catalog import read_cat
from .image import get_bgps_img


class Ds9FrameError(Exception):
    pass


class CoverageError(Exception):
    pass


class Inspector(object):
    rind_contour_color = 'yellow'
    bgps = read_cat('bgps_v210').set_index('v210cnum')
    v = 210
    inspect_cat = None

    def __init__(self, cnum):
        self.cnum = cnum
        self.set_coord()
        self.d = pyds9.DS9()
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
        #rind.verify(option='silentfix')
        rind[0].data[rind[0].data != self.cnum] = 0
        return rind

    def _get_flux(self):
        flux = get_bgps_img(self.cnum, exten='map20', v=210)
        #flux.verify(option='silentfix')
        return flux

    def set_coord(self):
        self.glonp = self.bgps.loc[self.cnum, 'glon_peak']
        self.glatp = self.bgps.loc[self.cnum, 'glat_peak']
        self.glonc = self.bgps.loc[self.cnum, 'glon_cen']
        self.glatc = self.bgps.loc[self.cnum, 'glat_cen']

    def show_source_contour(self):
        rind = self._get_rind()
        fid = self._get_next_frame_id()
        self.d.set('frame {0}'.format(fid))
        self.d.set_pyfits(rind)
        self.d.set('contour width 2')
        self.d.set('contour limits 0 1')
        self.d.set('contour nlevels 1')
        self.d.set('contour smooth 1')
        self.d.set('contour method block')
        self.d.set('contour')
        self.d.set('contour copy')
        self.d.set('frame 1')
        self.d.set('contour paste {0}'.format(self.rind_contour_color))
        self.d.set('frame delete {0}'.format(fid))

    def show_flux_contour(self, clevels, color='cyan', crop_coords=None):
        fid = self._get_next_frame_id()
        self.d.set('frame {0}'.format(fid))
        flux_img = self._get_flux()
        self.d.set_pyfits(flux_img)
        if crop_coords is not None:
            self.d.set('crop {0} wcs galactic degrees'.format(crop_coords))
        self.d.set('contour smooth 2')
        self.d.set('contour nvelels {0}'.format(len(clevels)))
        self.d.set('contour levels {' + ' '.join([str(i) for i in
                                                  clevels]) + '}')
        self.d.set('contour')
        self.d.set('contour copy')
        self.d.set('frame 1')
        self.d.set('contour paste {0}'.format(color))
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

    def cat_write(self, filen='inspect', flag=None):
        with open(filen + '.cat', 'a') as out_cat:
            if flag is None:
                flag = raw_input('flag {0:0>4d} : '.format(self.cnum))
            if flag == 'q':
                return
            out_cat.write('{0},{1}\n'.format(self.cnum, flag))


class HiGalInspector(Inspector):
    """
    Visually inspect HiGal cutouts in DS9.
    """
    # Directories and file names
    root_dir = '/home/svobodb/research/Data/HiGal/'
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
                 'der2y': 'der2y_source',
                 'der2y45': 'der2y45_source',
                 'source': 'source',
                 'mask_1': 'mask_1._source'}
    zlevel = 8
    max_scale = 1e4
    low_clevels = [0.1, 0.2, 0.4, 0.6, 0.8]
    high_clevels = [1.0, 2.0, 3.0, 4.0, 5.0]

    def __init__(self, cnum, img_type='source'):
        super(HiGalInspector, self).__init__(cnum)
        self.nmap = read_cat('hg70_index_map').set_index('v210cnum')
        self.cnum = cnum
        self.hcnum = self.nmap.loc[cnum, 'hg70_name']
        if img_type not in self.img_types.keys():
            raise ValueError('Invalid img_type: {0}.'.format(img_type))
        self.img_type = img_type
        self.filen = self._format_img_infile()

    def _format_img_infile(self):
        return self.img_file.format(cnum=self.hcnum, img_type=self.img_type,
                                    img_str=self.img_types[self.img_type])

    def _get_hg_img(self):
        self.img = fits.open(self.filen)
        self.img.verify(option='silentfix')
        self.img_min = self.img[0].data.min()
        self.img_max = self.img[0].data.max()
        self.xpix_max = self.img[0].data.shape[0]
        self.ypix_max = self.img[0].data.shape[1]
        self.d.set_pyfits(self.img)
        self.d.set('lock frame wcs')
        self.d.set('wcs sky galactic skyformat degrees')

    def set_scale(self):
        self.d.set('scale linear')
        if self.img_max > self.max_scale:
            self.d.set('scale mode zmax')
        else:
            self.d.set('scale mode zscale')

    def update_view(self, cnum, img_type='source'):
        self.cnum = cnum
        self.hcnum = self.nmap.loc[cnum, 'hg70_name']
        self.set_coord()
        self.img_type = img_type
        self.filen = self._format_img_infile()
        self.d.set('frame delete all')
        self.d.set('frame 1')
        self.view()

    def view(self):
        self._get_hg_img()
        self.zoom(self.zlevel)
        self.set_scale()
        crop_coords = '{0} {1} {2} {3}'.format(
            self.glonc, self.glatc, 0.2, 0.2)
        self.show_source_contour()
        self.show_flux_contour(clevels=self.low_clevels,
                               crop_coords=crop_coords)
        self.show_flux_contour(clevels=self.high_clevels, color='blue',
                               crop_coords=crop_coords)


class MipsgalInspector(Inspector):
    """
    Visually inspect HiGal cutouts in DS9.
    """
    # Directories and file names
    root_dir = '/home/svobodb/research/Data/MIPSGAL/'
    img_dir = os.path.join(root_dir, 'mosaics24')
    img_file = os.path.join(img_dir,
                            'MG{coord_str}_{img_type}.fits')
    # Properties
    img_types = {'covg': 'covg_024',
                 'mask': 'maskcube_024',
                 'err': 'std_024',
                 'source': '024'}
    zlevel = 8
    max_scale = 1e4
    low_clevels = [0.1, 0.2, 0.4, 0.6, 0.8]
    high_clevels = [1.0, 2.0, 3.0, 4.0, 5.0]

    def __init__(self, cnum, img_type='source'):
        super(MipsgalInspector, self).__init__(cnum)
        self.cnum = cnum
        if img_type not in self.img_types.keys():
            raise ValueError('Invalid img_type: {0}.'.format(img_type))
        self.img_type = img_type
        self.filen = self._format_img_infile()

    def _format_img_infile(self):
        glon = self.bgps.loc[self.cnum, 'glon_peak']
        glat = self.bgps.loc[self.cnum, 'glat_peak']
        psign = 'p'
        if glat < 0:
            psign = 'n'
        r_glon = int(round(glon))
        r_glat = int(round(glat))
        coord_str = '{glon:0>4d}{psign}{glat:0>3d}'.format(glon=r_glon,
                                                           psign=psign,
                                                           glat=r_glat)
        return self.img_file.format(coord_str,
                                    img_str=self.img_types[self.img_type])


def higal_group_inspect(filen='higal_inspect'):
    bgps = read_cat('bgps_v210').set_index('v210cnum')
    bgps = bgps.query('7.5 < glon_peak < 65')
    try:
        insp = pd.read_csv(filen + '.cat')
        cstart = insp.iloc[-1][0]
        bgps = bgps.loc[bgps.index > cstart]
    except IOError:
        pass
    # random cnum since has to start with one
    hgi = HiGalInspector(cnum=5000)
    for cnum in bgps.index:
        hgi.update_view(cnum=cnum)
        hgi.cat_write(filen=filen)


