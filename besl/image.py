"""
============
Image Module
============

Routines for images.

"""

import numpy as _np
from astropy.io import fits
from astropy import wcs
from .catalog import (read_bgps, select_bgps_field)
from .paths import all_paths as d


class Dirs(object):
    """
    Object to hold directories for interactive editing of paths.
    """
    def __init__(self):
        self.root_dir = d.root_dir
        self.bgps_dir = self.root_dir + 'BGPS/Images/{}/'
        self.bgps_img_filen = '{}_{}_13pca_{}.fits'
d = Dirs()


def get_bgps_img(identifier, exten, v=210):
    """
    Retrieve BGPS image file in astropy.io.fits instance. Only works for v2.

    Parameters
    ----------
    identifier : number
        BGPS catalog number of clump or a BGPS field string
    exten : string
        Image file extension name. Valid types:
        labelmask  -> source contours, label masks
        labelmap50 -> source contours, label masks for v1
        map20      -> default map
        map50      -> default map for v1
        medmap20   -> median map 20
        noisemap20 -> rms map
    v : number, default 2
        BGPS version number, valid [1, 2, 201, 'v2d', 210]. This only effects
        choice in label mask.

    Returns
    -------
    img : astropy.io.fits.HDUList
    """
    if v not in [101, 200, 201, '2d', 210]:
        raise ValueError('Invalid version: {0}'.format(v))
    if exten not in ['labelmask', 'labelmap50', 'map20', 'map50', 'medmap20',
                     'noisemap20']:
        raise ValueError('Incorrect exten: {}.'.format(exten))
    if v == '2d':
        exten = 'labelmask_deeper'
    ver_path = {101: 'v1.0.2', 200: 'v2.0.0', 201: 'v2.0.1', '2d': 'v2.0.1d',
                210: 'v2.1.0'}
    ver_init = {101: 'v1.0.2', 200: 'v2.0_ds2', 201: 'v2.0_ds2',
                '2d': 'v2.0_ds2', 210: 'v2.1_ds2'}
    cnum_col = {101: 'cnum', 200: 'cnum', 201: 'cnum',
                '2d': 'cnum', 210: 'v210cnum'}
    # cnum or field
    if isinstance(identifier, (float, int)):
        bgps = read_bgps(exten='none', v=v)
        c_index = _np.argwhere(bgps[cnum_col[v]] == identifier)[0][0]
        field = bgps.ix[c_index, 'field']
    elif isinstance(identifier, (str)):
        field = identifier
    else:
        raise ValueError('Improper identifier {0}.'.format(identifier))
    img = fits.open(d.bgps_dir.format(ver_path[v]) +
        d.bgps_img_filen.format(ver_init[v], field, exten))
    return img


def sample_img(img, coord):
    """
    Parameters
    ----------
    img : astropy.io.hdu.hdulist.HDUList
        Fits image.
    coord : tuple
        Coordinates in (lon, lat) form appropriate for native coordinates in
        the fits image header.

    Returns
    -------
    sample : number
        Returns `np.nan` if coordinate not present in image
    """
    img_wcs = wcs.WCS(img[0].header)
    # Convert coordinates to pixel values
    pix = _np.round(img_wcs.wcs_world2pix(coord[0], coord[1], 1)).astype(int)
    # Sample pixel value
    try:
        sample = img[0].data[pix[1], pix[0]]
    except:
        sample = _np.nan
    finally:
        return sample


def sample_bgps_img(lon, lat, exten='labelmask', v=210):
    """
    Retrieve a value from the BGPS images or labelmasks at a coordinate
    position.

    Parameters
    ----------
    lon : number
    lat : number
        Galactic coordinates in decimal degrees.

    Returns
    -------
    sample : number
        Returns `np.nan` if not contained in the BGPS bounds file.
    """
    # Get field identifier at coordinates
    field = select_bgps_field(lon=lon, lat=lat, coord_type='gal')
    if not isinstance(field, str):
        return _np.nan
    img = get_bgps_img(field, exten=exten, v=v)
    return sample_img(img=img, coord=(lon, lat))


class BgpsLib(object):
    """
    Container for BGPS images.
    """
    def __init__(self, exten='map20', v=210):
        """
        Parameters
        ----------
        exten : string, default 'map20'
            Image file extension name. Valid types:
            labelmask  -> source contours, label masks
            labelmap50 -> source contours, label masks for v1
            map20      -> default map
            map50      -> default map for v1
            medmap20   -> median map 20
            noisemap20 -> rms map
        v : number, default 210
            BGPS version number
        """
        self.exten = exten
        self.v = v
        self.bgps = read_bgps()
        self.fields = self.bgps.field.unique()
        self._images = {}

    def read_images(self):
        for field in self.fields:
            img = get_bgps_img(field, exten=self.exten, v=self.v)
            self._images[field] = img

    def get_hdu(self, field):
        return self._images[field]

    def get_img(self, field):
        return self._images[field][0].data

    def get_hdr(self, field):
        return self._images[field][0].header


