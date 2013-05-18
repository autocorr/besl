"""
============
Image Module
============

Routines for images.

"""

import pyfits
import catalog
import numpy as _np

class Dirs(object):
    """
    Object to hold directories for interactive editing of paths.
    """
    def __init__(self):
        self.root_dir = '/mnt/eld_data/'
        self.bgps_dir = self.root_dir + 'BGPS/Images/v2.0.0/'
        self.bgps_img_filen = 'v2.0_ds2_{}_13pca_{}.fits'
d = Dirs()

def get_bgps_img(cnum, exten):
    """
    Retrieve BGPS image file in pyfits format. Only works for v2.

    Parameters
    ----------
    cnum : number
        BGPS catalog number of clump
    exten : string
        Image file extension name. Valid types:
        labelmask  -> source contours, label masks
        map20      -> default map
        medmap20   -> median map 20
        noisemap20 -> rms map

    Returns
    -------
    img : pyfits.HDUList
    """
    if exten not in ['labelmask', 'map20', 'medmap20', 'noisemap20']:
        raise ValueError('Incorrect exten: {}.'.format(exten))
    bgps = catalog.read_bgps(exten='none')
    c_index = _np.argwhere(bgps.cnum == cnum)[0][0]
    field = bgps.ix[c_index, 'field']
    img = pyfits.open(d.bgps_dir + d.bgps_img_filen.format(field, exten))
    return img


