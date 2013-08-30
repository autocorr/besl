"""
==================
PPV Grouping Plots
==================

Plotting routines for the PPV Grouping and velocity visualization.

"""

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from ..catalog import read_bgps_vel, read_bgps_bounds


def velocity_color_scatter(field=None, coords=None):
    """
    Parameters
    ----------
    field : str, default None
        If left None, then specify longitude range
    coords : list

    Returns
    -------
    (fig, ax) : tuple
    """
    # data
    bounds = read_bgps_bounds()
    bgpv = read_bgps_vel()
    # limits
    ii = bounds[bounds.field == field].index[0]
    glon_min = bounds.ix[ii, 'glon_min']
    glon_max = bounds.ix[ii, 'glon_max']
    glat_min = bounds.ix[ii, 'glat_min']
    glat_max = bounds.ix[ii, 'glat_max']
    mask = (bgpv.glon_peak > glon_min) & (bgpv.glon_peak < glon_max) & \
           (bgpv.glat_peak > glat_min) & (bgpv.glat_peak < glat_max) & \
           (bgpv.vlsr_f > 0)
    glons = bgpv[mask]['glon_peak'].values
    glats = bgpv[mask]['glat_peak'].values
    velos = bgpv[mask]['vlsr'].values
    # figure
    fig, ax = plt.subplots()
    ax.grid()
    sc = ax.scatter(glons, glats, c=velos, cmap=cm.jet, vmin=velos.min(),
                    vmax=velos.max())
    cb = plt.colorbar(sc, ax=ax)
    # limits
    ax.set_xlim([glon_min, glon_max])
    ax.set_ylim([glat_min, glat_max])
    ax.axis('equal')
    plt.draw()
    return fig, ax


