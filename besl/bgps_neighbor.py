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
    return

