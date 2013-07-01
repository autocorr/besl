"""
==========================
Kinetic Temperatures Plots
==========================

Plots for kinetic temperatures (T_K) and other line properties derived from
NH3 observations.

"""
# TODO add plots for line widths and intensities

import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from .. import catalog
from .. import util
import ipdb as pdb

def ktemp_hist(bgps=[]):
    """
    Plot a histogram of all kinetic temperatures
    """
    bgps = util.bgps_import_check(bgps)
    # Select data
    bgps['nh3_snr11'] = bgps['nh3_pk11'] / bgps['nh3_noise11']
    bgps['nh3_snr22'] = bgps['nh3_pk22'] / bgps['nh3_noise22']
    bgps['nh3_snr33'] = bgps['nh3_pk33'] / bgps['nh3_noise33']
    ktemps = bgps[(bgps.nh3_snr11 > 4) & (bgps.nh3_snr22 > 4)]['nh3_tkin'].values
    ktemps = ktemps[(ktemps > 5) & (ktemps < 150)]
    # Plot settings
    hist_kwargs = {'bins': _np.logspace(_np.log10(ktemps.min()),
        _np.log10(ktemps.max()), 40), 'histtype': 'stepfilled', 'color': 'gray',
        'log': True}
    fig = _plt.figure(figsize=(6,4))
    fig.subplots_adjust(bottom=0.12)
    ax = fig.add_subplot(111)
    # Begin plot
    ax.hist(ktemps, **hist_kwargs)
    #ax.plot([40, 40], [0.8, 100], 'k-.')  # 40 K for saturate 11/22
    ax.plot(_np.median(ktemps) * _np.ones(2), [0.8, 100], 'k--')
    ax.set_xlabel(r'$T_{\rm K} \ \ [{\rm K}]$')
    ax.set_ylabel(r'$N$')
    ax.set_xscale('log')
    ax.set_ylim(bottom=0.8)
    return ax


