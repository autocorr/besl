"""
==================
PPV Grouping Plots
==================

Plotting routines for the PPV Grouping and velocity visualization.

"""

import cPickle as pickle
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.colors import LogNorm
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

def tree_params(filen='obj_props', out_filen='ppv_grid', log_Z=False):
    """
    Read in the pickled tree parameter dictionary and plot the containing
    parameters.

    Parameters
    ----------
    filen : str
        File name of pickled reduced property dictionary.
    out_filen : str
        Basename of plots, the key of the object dictionary is appended to the
        filename.
    log_Z : bool
        Create plots with logarithmic Z axis
    """
    obj_dict = pickle.load(open(filen + '.pickle', 'rb'))
    X = obj_dict['velo']
    Y = obj_dict['angle']
    obj_dict['reward'] = np.log10(obj_dict['new_kdar_assoc']) / obj_dict['conflict_frac']
    params = [(k, v) for k, v in obj_dict.iteritems()
              if k not in ['velo', 'angle']]
    clevels = [0.05, 0.10, 0.15]
    for key, Z in params:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if log_Z:
            Z = np.log10(Z)
            key += '_(log)'
        pc = ax.pcolor(X, Y, Z, cmap=cm.cubehelix, vmin=Z.min(), vmax=Z.max())
        cb = plt.colorbar(pc, ax=ax)
        # Contours for conflict frac
        cn = ax.contour(X, Y, obj_dict['conflict_frac'], levels=clevels,
                        colors='k', linewidth=2)
        plt.setp(cn.collections,
                 path_effects=[PathEffects.withStroke(linewidth=2,
                 foreground='w')])
        cl = ax.clabel(cn, fmt='%1.2f', inline=1, fontsize=12,
                       use_clabeltext=True)
        plt.setp(cl, path_effects=[PathEffects.withStroke(linewidth=2,
                 foreground='w')])
        # Labels
        ax.set_xlabel(r'$v \ \ [{\rm km \ s^{-1}}]$')
        ax.set_ylabel(r'$\theta \ \ [^{\circ}]$')
        ax.set_title(' '.join(key.split('_')))  # latex doesn't like _
        # Limits
        ax.set_xlim([X.min(), X.max()])
        ax.set_ylim([Y.min(), Y.max()])
        # Save
        plt.savefig(out_filen + '_' + key + '.pdf')
        plt.savefig(out_filen + '_' + key + '.png')

def conflict_hist(filen='good_cluster', out_filen='conflict_hist'):
    """
    Read in the pickled `besl.ppv_group.ClusterDBSCAN` instance and plot a
    histogram of the number of clusters with a given number of nodes.

    Parameters
    ----------
    filen : str
        File name of pickled reduced property dictionary.
    out_filen : str
        Basename of plots, the key of the object dictionary is appended to the
        filename.

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.AxesSubplot
    """
    # Read in data
    obj = pickle.load(open(filen + '.pickle', 'rb'))
    conflicts, agrees, singles = [], [], []
    for cid, params in obj.cluster_nodes.iteritems():
        nodes = params[0]
        conflict_flag = params[2]
        if conflict_flag == 1:
            agrees.append(len(nodes))
        elif (conflict_flag == 2) & (cid != -1):
            conflicts.append(len(nodes))
        elif (conflict_flag == 0):
            singles.append(len(nodes))
    nmax = max([max(c) for c in [conflicts, agrees, singles]])
    bins = np.arange(0, nmax + 1)
    kwargs = {'bins': bins, 'log': True}
    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(singles, color='grey', linewidth=2, histtype='stepfilled',
            **kwargs)
    ax.hist(conflicts, color='red', hatch='/', histtype='step', **kwargs)
    ax.hist(agrees, color='green', hatch='\\', histtype='step', **kwargs)
    ax.set_xlabel(r'${\rm Nodes}$')
    ax.set_ylabel(r'$N$')
    ax.set_yscale('log')
    ax.set_ylim(bottom=0.8)
    # Save
    plt.savefig(out_filen + '.pdf')
    plt.savefig(out_filen + '.png')
    return fig, ax


