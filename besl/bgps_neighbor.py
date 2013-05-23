"""
=====================
BGPS Nearest Neighbor
=====================

Match and analyze BGPS nearest neighbors in label masks.

"""

import pyfits
import numpy as _np
import pandas as _pd
import catalog, image

def find_clump_neighbors(cnum, v=201):
    """
    Determine the neighbors of a clump by if there is adjacent overlap in the
    label masks.

    Parameters
    ----------
    cnum : number
        BGPS source catalog number
    v : number, default 201
        BGPS version number. Available options:
            2   -> BGPS v2.0.0
            201 -> BGPS v2.0.1

    Returns
    -------
    neighbors : list
        List of neighbor's catalog numbers. If no neighbors, returns an empty
        list.
    """
    if v not in [2, 201]:
        raise ValueError
    rind = image.get_bgps_img(cnum, exten='labelmask', v=v)
    # find perimeter pixels
    perim_pixels = _np.argwhere(
         (rind[0].data[1:-1,1:-1] == cnum) &
        ((rind[0].data[0:-2,1:-1] != cnum) |
         (rind[0].data[2:  ,1:-1] != cnum) |
         (rind[0].data[1:-1,0:-2] != cnum) |
         (rind[0].data[1:-1,2:  ] != cnum)))
    # look up, down, left, right
    udlr_pixels = []
    for i, j in [(0,1), (0,-1), (1,0), (-1,0)]:
        udlr_pixels.append(
            _np.dstack([perim_pixels[:,0] + i, perim_pixels[:,1] + j])) # zip
    # re-sample label masks
    neighbors = []
    for pixels in udlr_pixels:
        neighbors.append(rind[0].data[pixels])
    return list(_np.unique(neighbors))

