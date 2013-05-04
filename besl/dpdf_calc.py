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

import os as _os
import numpy as _np
import pandas as _pd
import matplotlib.pyplot as _plt
from besl.catalog import read_dpdf, read_emaf_dist, read_oh94_dust
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
    xsamples : np.array
    ysamples : np.array
        Arrays of sampled x and evaluated y values
    """
    if len(lims) != 4:
        raise ValueError('lims must be length 4.')
    if len(x) != len(y):
        raise ValueError('x and y must be equal in length')
    xmin, xmax, ymin, ymax = lims
    N = len(x)
    fn = interp1d(x, y, kind='linear')
    xsamples = []
    # sample distribution
    while len(xsamples) < nsample:
        rands = _np.random.uniform(low=(xmin, ymin), high=(xmax, ymax),
            size=(nsample, 2))
        xsamples = _np.append(xsamples, rands[rands[:,1] < fn(rands[:,0]), 0])
    # clip and evaluate samples
    xsamples = xsamples[0:nsample]
    ysamples = fn(xsamples)
    return xsamples, ysamples

def clump_dust_mass(dist, snu=1, tkin=20., nu=2.725e11):
    """
    Calculate the dust mass from the distance, specific flux density, and
    kinetic temperature.

    Parameters
    ----------
    dist : array-like
        Distance in pc
    snu : number, default 1
        Specific flux density at `nu' in Jy
    tkin : number, default
        Kinetic Temperature in K
    nu : number, defualt
        Frequency in GHz

    Returns
    -------
    mdust : np.array
    """
    # TODO add blackbody function
    from besl.units import cgs
    bnu = planck_fn(nu, tkin, freq=True)
    oh5 = read_oh94_dust(model_type='thick', modeln=0)
    kapp = oh5(2.)
    return (snu * dist**2) / (kapp * bnu)

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
    x = _np.arange(1000) * 20. + 20.
    y = dpdf[5].data[0]
    lims = [x.min(), x.max(), 0, 1]
    samples = mc_sampler_2d(x, y, lims=lims, nsample=1e4)
    fig = _plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(samples, bins=_np.linspace(lims[0], lims[1], 200), linewidth=2,
        histtype='stepfilled', normed=True, alpha=0.5, color='black')
    ax.plot(x, y / 20, 'k-', linewidth=2)
    ax.set_xlabel(r'$D \ \ [{\rm pc}]$')
    ax.set_ylabel(r'${\rm Probability}$')
    ax.set_ylim([0, y.max() * 1.1 / 20.])
    _plt.savefig('dpdf_test_sampling.pdf', format='pdf')
    return ax

def print_dpdf_outfiles(out_dir='dpdf_ascii', v=2):
    """
    Print out ascii files for the posterior DPDF.

    Parameters
    ----------
    out_dir : string, default 'dpdf_ascii'
        Directory to place dpdf files
    v : number, default 2
        BGPS version number
    """
    if not os.path.exists(out_dir):
        raise Exception('Out path does not exist.')
    if v not in [1, 2]:
        raise ValueError('Incorrect version.')
    dpdf = read_dpdf(v=v)
    flags = _pd.DataFrame(dpdf[1].data)
    for i, row in enumerate(dpdf[5].data):
        _np.savetxt(out_dir + '/v{0}_{1:0>4d}.txt'.format(v,
            flags.CNUM.iloc[i]), row, fmt='%.10e')
    return
