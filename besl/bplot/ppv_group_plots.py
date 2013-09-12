"""
==================
PPV Grouping Plots
==================

Plotting routines for the PPV Grouping and velocity visualization.

"""

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ..catalog import read_bgps_vel


def velocity_color_scatter(coords):
    """
    Plot a scatter plot of velocity positions with velocities represented by
    differently colored points.

    Parameters
    ----------
    coords : list
        Coordinate list of [glon_min, glon_max, glat_min, glat_max] in decimal
        degrees.

    Returns
    -------
    (fig, ax) : tuple
    """
    assert len(coords) == 4
    # data
    bgpv = read_bgps_vel()
    # limits
    glon_min, glon_max, glat_min, glat_max = coords
    mask = (bgpv.glon_peak > glon_min) & (bgpv.glon_peak < glon_max) & \
           (bgpv.glat_peak > glat_min) & (bgpv.glat_peak < glat_max) & \
           (bgpv.vlsr_f > 0)
    glons = bgpv[mask]['glon_peak'].values
    glats = bgpv[mask]['glat_peak'].values
    velos = bgpv[mask]['vlsr'].values
    # figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    sc = ax.scatter(glons, velos, glats, c=velos, cmap=cm.jet, vmin=velos.min(),
               vmax=velos.max())
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
    # limits
    ax.set_xlim([glon_min, glon_max])
    ax.set_ylim([velos.min(), velos.max()])
    ax.set_zlim([glat_min, glat_max])
    ax.invert_yaxis()
    # labels
    ax.set_xlabel(r'$\ell \ \ [^{\circ}]$')
    ax.set_ylabel(r'$v \ \ [{\rm km \ s^{-1}}]$')
    ax.set_zlabel(r'$b \ \ [^{\circ}]$')
    plt.draw()
    return fig, ax

def tree_params():
    """
    """
    pass


