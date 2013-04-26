"""
================
DPDF Monte Carlo
================

Routines for handling EMAF DPDF's. Citation: Ellsworth-Bowers et al. (2013).
"""
# TODO read DPDFs
# create functions for:
#   dust masses
#   sizes / radius
#   bolometric luminosity
#   isotropic water maser luminosity
# monte carlo sampler

import pyfits
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from besl.catalog import read_dpdf, read_emaf_dist
from scipy.interpolate import interp1d

def mc_sampler_2d(x, y, lims=[0,1,0,1], nsample=1e3):
    """
    Monte Carlo sampler.

    Parameters
    ----------
    x, y : np.array
        Function to sample, must be same shape. Interpolated using
        scipy.interpolate.interp1d.
    lims : list, default [0,1,0,1]
        Limits to uniformly sample from [xmin, xmax, ymin, ymax]
    nsample : number, default 1e3
        Number of samples to return

    Returns
    -------
    samples : np.array
        Array of sampled values
    """
    # TODO put in mathf module
    if len(lims) != 4:
        raise ValueError('lims must be length 4.')
    if len(x) != len(y):
        raise ValueError('x and y must be same length.')
    xmin, xmax, ymin, ymax = lims
    N = len(x)
    fn = interp1d(x, y, kind='linear')
    samples = []
    while len(samples) < nsample:
        rands = np.random.uniform(low=(xmin, ymin), high=(xmax, ymax),
            size=(nsample, 2))
        samples = np.append(samples, rands[rands[:,1] < fn(rands[:,0]), 0])
    return samples[0:nsample]

def sample_distrib(x, y, samples):
    """
    Sample a distribution at selected positions.

    Parameters
    ----------
    x, y : np.array
        Data to sample, will

    Returns
    -------

    """
    # TODO put in mathf module
    pass

def clump_dust_mass():
    pass

def clump_radius():
    pass

def clump_dust_luminosity():
    pass

def clump_line_luminosity():
    pass

def plot_dpdf_sampling(n=200):
    """
    Plot a Monte Carlo sampling of a DPDF
    """
    dpdf = read_dpdf()
    x = np.arange(1000) * 20. + 20.
    y = dpdf[5].data[0]
    lims = [x.min(), x.max(), 0, 1]
    samples = mc_sampler_2d(x, y, lims=lims, nsample=1e4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(samples, bins=np.linspace(lims[0], lims[1], 200), linewidth=2,
        histtype='stepfilled', normed=True, alpha=0.5, color='black')
    ax.plot(x, y / 20, 'k-', linewidth=2)
    ax.set_xlabel(r'$D \ \ [{\rm pc}]$')
    ax.set_ylabel(r'${\rm Probability}$')
    ax.set_ylim([0, y.max() * 1.1 / 20.])
    plt.show()
    return ax
