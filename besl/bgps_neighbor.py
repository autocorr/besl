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
    neighbors : np.array
        Array of unique neighbor's catalog numbers. If no neighbors, returns an
        empty array.
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
    neighbors = _np.unique(neighbors)
    neighbors = neighbors[neighbors != 0]
    if _np.any(neighbors < 0):
        raise ValueError('Negative cnum found')
    return neighbors

def broadcast_kdar():
    # TODO
    # for clumps with a dpdf
    #   get list of neighbor clumps
    #   visit each neighbor, check if it has a velocity, and if that velocity
    #     is close enough to the DPDF clumps velocity
    #   update visited clump list
    #   once all adjacent neighbors have been visited, look at the neighbors of
    #   the clumps in the visited list but only if that neighbor is not already
    #   in the visited list
    #   once out of clumps to visit
    pass

def num_of_neighbors(v=201):
    """
    Calculate the number of neighbors a clump has based on the label masks.

    Parameters
    ----------
    v : number, default 201
        BGPS version number. Available options:
            2   -> BGPS v2.0.0
            201 -> BGPS v2.0.1

    Returns
    -------
    bgps : pd.DataFrame
        dataframe with 'neighbors_n' column
    """
    bgps.catalog.read_bgps(exten='none', v=v)
    bgps['neighbors_n'] = _np.nan
    for i in xrange(bgps.shape[0]):
        cnum, field = bgps.ix[i, ['cnum', 'field']]
        neighbors = find_clump_neighbors(cnum, v=v)
        bgps.ix[i, 'neighbors_n'] = len(neighbors)
    bgps = bgps[['name', 'neighbors_n']]
    return bgps


