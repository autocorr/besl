"""
===============
BGPS Clump Area
===============

Calculate the area of a BGPS clump on the sky.

"""

import catalog, mathf, image
import numpy as _np
import pandas as _pd

def bgps_rind_area(cnum):
    """
    Calculate the area of a BGPS clump from the number of pixels in the
    rind map / labelmask.

    Parameters
    ----------
    cnum : number
        BGPS v2 catalog number

    Returns
    -------
    area : number
        Area of clump in square arcseconds
    """
    bgps = catalog.read_bgps()
    field_cnum = bgps[bgps.cnum == cnum]['field_cnum'].values[0]
    rind = image.get_bgps_img(cnum, exten='labelmask')
    clump_pixels = _np.argwhere(rind[0].data == field_cnum)
    area = clump_pixels.shape[0] * 7.5**2
    return area

def bgps_rind_perim_pix(cnum):
    """
    Calculate the number of perimeter pixels for a BGPS clump from the rind map
    / labelmask.

    Parameters
    ----------
    cnum : number
        BGPS v2 catalog number

    Returns
    -------
    perim : number
        Perimeter of ellipse in square arcseconds
    """
    bgps = catalog.read_bgps()
    field_cnum = bgps[bgps.cnum == cnum]['field_cnum'].values[0]
    rind = image.get_bgps_img(cnum, exten='labelmask')
    perim_pixels = _np.argwhere(
         (rind[0].data[1:-1,1:-1] == cnum) &
        ((rind[0].data[0:-2,1:-1] != cnum) |
         (rind[0].data[2:  ,1:-1] != cnum) |
         (rind[0].data[1:-1,0:-2] != cnum) |
         (rind[0].data[1:-1,2:  ] != cnum)))
    return perim_pixels.shape[0]


