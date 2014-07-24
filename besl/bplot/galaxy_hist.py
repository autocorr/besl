"""
===================
Galaxy 2D Histogram
===================

Rourtines for plotting a 2D rectangular binned histrogram over a large range of
Galactic longitude using coarse binning.

"""

import copy
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import patheffects
from ..dpdf_calc import evo_stages


def bgps_stages(bin_int=2, bgps=None, verbose=False):
    """
    bin_int : number
        Bin-size interval, number of sub-divisions per degree. Square bins are
        used with equal length in the xy plane.
    """
    stages, anno_labels = evo_stages(bgps=bgps)
    xlims = [10, 60]
    ylims = [-1, 1]
    xbins = np.linspace(xlims[0], xlims[1], 60 * bin_int + 1)
    ybins = np.linspace(ylims[0], ylims[1], 2 * bin_int + 1)
    stages_labels = [i.values()[0] for i in anno_labels]
    # color scale
    cmap = copy.copy(cm.summer)
    cmap._init()
    cmap.set_bad('w', 0)
    cmap._lut[0,:] = 0
    # create plot
    fig, axes = plt.subplots(nrows=len(stages), ncols=1, sharex=True)
    for i, ax in enumerate(axes):
        df = stages[i]
        lons = df.glon_peak.values
        lats = df.glat_peak.values
        # make histogram
        Hist, xedges, yedges = np.histogram2d(lons, lats, bins=[xbins, ybins])
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        ax.fill([xedges[0], xedges[-1], xedges[-1], xedges[0]],
                [yedges[0], yedges[0], yedges[-1], yedges[-1]], fill=False,
                hatch='//', zorder=1, color='0.5', linewidth=1.5,
                antialiased=True)
        im = ax.imshow(Hist.T, extent=extent, aspect='auto', cmap=cmap,
                       vmax=Hist.max(), interpolation='nearest', zorder=2)
        ax.invert_xaxis()
        # annotate
        name_txt = ax.annotate(stages_labels[i], xy=(0.025, 0.65),
                               xycoords='axes fraction', fontsize=12,
                               weight='bold', zorder=3)
        name_txt.set_path_effects([patheffects.withStroke(linewidth=2,
            foreground='w')])
        stage_counts = df.shape[0]
        count_txt = ax.annotate(stage_counts, xy=(0.025, 0.15),
                                xycoords='axes fraction', fontsize=12,
                                weight='bold', zorder=3)
        count_txt.set_path_effects([patheffects.withStroke(linewidth=2,
            foreground='w')])
        # ticks
        ax.set_yticks([-0.5, 0, 0.5])
        ax.minorticks_on()
        ax.tick_params(axis='y', which='minor', left='off', right='off')
        # colorbar
        #cb = plt.colorbar(im, ax=ax)
        if verbose:
            print 'Max cell in histogram: {0}'.format(Hist.max())
    # plot axes labels
    axes[-1].set_xlabel(r'$\ell \ \ [^{\circ}]$')
    axes[-1].set_ylabel(r'$b \ \ [^{\circ}]$')
    # save
    plt.subplots_adjust(hspace=0.1)
    plt.savefig('bgps_galaxy_hist.pdf')
    print '-- bgps_galaxy_hist.pdf written'
    return fig, axes
