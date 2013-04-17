"""
Authors: Brian Svoboda
Date: 2013, 4, 11

Plotting routines for BGPS observation boundaries
"""

import numpy as _np
import matplotlib.pyplot as _plt

def aitoff_proj():
    # TODO broken, use newer verson of pywcsgrid2 to make projection
    """
    Plot BGPS boundaries in an Aitoff projection.

    Returns
    -------
    fig : matplotlib figure object
    ax : matplotlib axes object
    """
    from besl.catalog import read_bgps_bounds
    import pywcsgrid2
    # read in catalog
    bgpsb = read_bgps_bounds()
    fig = _plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='aitoff')
    ax.invert_xaxis()
    ax.set_xlim([0, 360])
    ax.grid()
    lon_max_mask = bgpsb.glon_max > 180
    lon_min_mask = bgpsb.glon_max > 180
    bgpsb.glon_max[lon_max_mask] -= 360
    bgpsb.glon_min[lon_min_mask] -= 360
    bgpsb.glon_max = _np.deg2rad(bgpsb.glon_max)
    bgpsb.glon_min = _np.deg2rad(bgpsb.glon_min)
    bgpsb.glat_max = _np.deg2rad(bgpsb.glat_max)
    bgpsb.glat_min = _np.deg2rad(bgpsb.glat_min)
    ax.plot(bgpsb['glon_max'], bgpsb['glat_max'], 'r.')
    ax.plot(bgpsb['glon_min'], bgpsb['glat_min'], 'b.')
    return [fig, ax, bgpsb]
