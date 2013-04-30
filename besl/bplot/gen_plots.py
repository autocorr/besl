"""
=============
General Plots
=============

Plotting routines for non-specialized applications.

"""

import matplotlib.pyplot as plt
from bconfig import *

def bhexbin(x, y, labels=[], lims=[], out_filen=None):
    """
    Plot a single 2-dimensional histogram with hexagon binning.

    Parameters
    ----------
    x, y : array-like
    labels : list
        List of x axis, y axis, and colorbar labels
    lims : list
        List of x min, x max, y min, y max plot limits
    out_filen : string

    Returns
    -------
    ax : matplotlib.Axes
    """
    if len(x) != len(y):
        raise ValueError
    if (type(labels) != list) & (type(lims) != list):
        raise TypeError
    if len(labels) == 0:
        xlabel, ylabel, clabel = 3 * ['']
    else:
        xlabel, ylabel, clabel = labels
    if len(lims) == 0:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        lims = [xmin, xmax, ymin, ymax]
    fig = plt.figure(figsize=(6,4.5))
    ax = fig.add_subplot(111)
    hx = ax.hexbin(x, y, **bhexbin_kwargs)
    cb = plt.colorbar(hx, ax=ax)
    cb.set_label(r'$ {}$'.format(clabel))
    ax.set_xlabel(r'$ {}$'.format(xlabel))
    ax.set_ylabel(r'$ {}$'.format(ylabel))
    ax.set_xlim([lims[0], lims[1]])
    ax.set_ylim([lims[2], lims[3]])
    if out_filen is not None:
        plt.savefig(out_filen + '.png', dpi=300)
    return ax

def bhist():
    pass
