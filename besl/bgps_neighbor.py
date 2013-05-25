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
import ipdb as pdb

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
         (rind[0].data[1:-1,2:  ] != cnum))) + 1 # for pyfits index
    # look up, down, right, left
    udlr_pixels = []
    for i, j in [(0,1), (0,-1), (1,0), (-1,0)]:
        udlr_pixels.append(
            _np.array([perim_pixels[:,0] + i, perim_pixels[:,1] + j]).T) # zip
    # re-sample label masks
    neighbors = []
    for pixels in udlr_pixels:
        neighbors.append(rind[0].data[pixels[:,0], pixels[:,1]])
    neighbors = _np.unique(_np.ravel(neighbors))
    neighbors = neighbors[(neighbors != 0) & (neighbors != cnum)]
    if _np.any(neighbors < 0):
        raise ValueError('Negative cnum found')
    return neighbors

def select_good_neighbors(bgps, cnum, hco_v, visited):
    """
    Select neighbors of clump that satisfy the criteria:
      * Have HCO+ flag 1 or 3
      * HCO+ velocity closer than 3.5 km / s
      * Do not have a KDAR 'N', 'F', or 'T'
      * Have not already been visited.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS all catalog
    cnum : number
        Current clump v2.0.1 catalog number
    visited : array-like
        List of visited v2.0.1 catalog numbers

    Returns
    -------
    good_neighbors : list
    """
    all_neighbors = find_clump_neighbors(cnum)
    good_neighbors = bgps.ix[
        (bgps.v201cnum.isin(all_neighbors)) &
        (_np.logical_not(bgps.v201cnum.isin(visited))) &
        (_np.logical_not(bgps.KDAR.isin(['N','F','T']))) &
        (_np.abs(bgps.hco_v - hco_v) < 3.5) &
        (bgps.hco_f.isin([1,3])), 'v201cnum'].values
    return list(good_neighbors)

def broadcast_kdar(bgps=[], verbose=False):
    """
    Calculate the nearest neighbors to the clumps with KDAR resolved by EMAF
    in Ellsworth-Bowers et al. (2013). Only v2.0.1 names supported.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS catalog dataframe
    verbose : bool, default False

    Returns
    -------
    bgps : pd.DataFrame
        BGPS catalog with added columns:
            neighbor_KDAR : KDAR from a spanning DPDF clump
            neighbor_dML  : dML from a spanning DPDF clump
    """
    # TODO
    # - mark special flag or raise an exception if conflicting KDAR within a
    #   complex
    # - call special routine if a found neighbor with a KDAR is found that will
    #   step each DPDF clump in a complex in parallel order
    if type(verbose) is not bool:
        raise TypeError
    if len(bgps) == 0:
        bgps = catalog.read_bgps(exten='all')
    bgps['neighbor_KDAR'] = 'null'
    bgps['neighbor_dML'] = _np.nan
    # visit DPDF clumps
    for i in bgps[bgps.KDAR.isin(['T','N','F'])].index:
        # current DPDF clump properties
        dpdf_cnum = bgps.ix[i, 'v201cnum']
        kdar = bgps.ix[i, 'KDAR']
        dML = bgps.ix[i, 'dML']
        hco_v = bgps.ix[i, 'hco_v']
        visited = [dpdf_cnum]
        neighbors = select_good_neighbors(bgps, dpdf_cnum, hco_v, visited)
        if verbose:
            print '\n-- DPDF clump : {}'.format(dpdf_cnum)
        for neighbor_cnum in neighbors:
            visited.append(neighbor_cnum)
            # update flags for current clump
            bgps['neighbor_KDAR'][bgps.v201cnum == neighbor_cnum] = kdar
            bgps['neighbor_dML'][bgps.v201cnum == neighbor_cnum] = dML
            # check for new clumps
            new_neighbors = select_good_neighbors(bgps, neighbor_cnum,
                hco_v, visited)
            neighbors.extend(new_neighbors)
            if verbose:
                print '.',
    bgps['neighbor_KDAR'][bgps.neighbor_KDAR == 'null'] = _np.nan
    return bgps

def num_of_neighbors(v=201, verbose=False):
    """
    Calculate the number of neighbors a clump has based on the label masks.

    Parameters
    ----------
    v : number, default 201
        BGPS version number. Available options:
            2   -> BGPS v2.0.0
            201 -> BGPS v2.0.1
    verbose : Bool, default False

    Returns
    -------
    bgps : pd.DataFrame
        dataframe with 'neighbors_n' column
    """
    bgps = catalog.read_bgps(exten='none', v=v)
    bgps['neighbors_n'] = _np.nan
    for i in xrange(bgps.shape[0]):
        cnum, field = bgps.ix[i, ['cnum', 'field']]
        neighbors = find_clump_neighbors(cnum, v=v)
        bgps.ix[i, 'neighbors_n'] = len(neighbors)
        if verbose:
            print '-- {0:>4d} : {1:>4d}'.format(cnum, len(neighbors))
    bgps = bgps[['name', 'neighbors_n']]
    return bgps


