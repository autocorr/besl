"""
============
Image Module
============

Routines for images.

"""

import numpy as _np
from .catalog import (read_bgps, select_bgps_field)
from astropy.io import fits
from astropy import wcs


class Dirs(object):
    """
    Object to hold directories for interactive editing of paths.
    """
    def __init__(self):
        self.root_dir = '/mnt/eld_data/'
        self.bgps_dir = self.root_dir + 'BGPS/Images/{}/'
        self.bgps_img_filen = '{}_{}_13pca_{}.fits'
d = Dirs()


def get_bgps_img(identifier, exten, v=201):
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
    if v not in [1, 2, 201, '2d', 210]:
        raise ValueError('Invalid version: {0}'.format(v))
    if exten not in ['labelmask', 'labelmap50', 'map20', 'map50', 'medmap20',
                     'noisemap20']:
        raise ValueError('Incorrect exten: {}.'.format(exten))
    if v == '2d':
        exten = 'labelmask_deeper'
    ver_path = {1: 'v1.0.2', 2: 'v2.0.0', 201: 'v2.0.1', '2d': 'v2.0.1d',
                210: 'v2.1.0'}
    ver_init = {1: 'v1.0.2', 2: 'v2.0_ds2', 201: 'v2.0_ds2',
                '2d': 'v2.0_ds2', 210: 'v2.0_ds2'}
    # cnum or field
    if isinstance(identifier, (float, int)):
        bgps = read_bgps(exten='none', v=v)
        c_index = _np.argwhere(bgps.cnum == identifier)[0][0]
        field = bgps.ix[c_index, 'field']
    elif isinstance(identifier, (str)):
        field = identifier
    else:
        raise ValueError('Improper identifier {0}.'.format(identifier))
    img = fits.open(d.bgps_dir.format(ver_path[v]) +
        d.bgps_img_filen.format(ver_init[v], field, exten))
    return img


def sample_bgps_img(coord):
    """
    Retrieve a value from the BGPS images or labelmasks at a coordinate
    position.

    Parameters
    ----------
    coord : tuple
        Galactic coordinates in decimal degrees. Formed in (lon, lat).

    Returns
    -------
    sample : number
        Returns `None` if not contained in the BGPS.
    """
    # Get field of image
    # If no field, return None
    # Use field to query `get_bgps_image` and get image
    # get wcs information
    # convert coordinates to pixel values
    # sample pixel value and return
    pass
