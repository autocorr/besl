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
import catalog
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
    # TODO add planck function
    from besl.units import cgs
    #bnu = planck_fn(nu, tkin, freq=True)
    oh5 = catalog.read_oh94_dust(model_type='thick', modeln=0)
    kapp = oh5(2.)
    return (snu * dist**2) / (kapp * bnu)

def clump_simple_dust_mass(dist, snu11, tkin=20.):
    """
    Calculate the dust mass in Solar masses based on an OH5 dust opacity model
    at 1.1mm, as used by the BGPS.

    Parameters
    ----------
    dist : array-like
        Distance in pc
    snu11 : number
        1.1 mm flux density in Jy
    tkin : number, default
        Kinetic Temperature in K

    Returns
    -------
    mdust_mean : number
    mdust_err : number
    """
    mdust_dist = 14.26 * (_np.exp(13.01 / tkin) - 1) * snu11 * (dist * 1e-3)**2
    mdust_mean = _np.mean(mdust_dist)
    mdust_err = _np.std(mdust_dist)
    return mdust_mean, mdust_err

def clump_surface_area(dist, maj_axis, min_axis):
    """
    Calculate the clump surface area in square pc based on the BGPS ellipse.

    Parameters
    ----------
    dist : array-like
        Distance in pc
    maj_axis : number
        Semi-major axis in arcseconds
    min_axis : number
        Semi-minor axis in arcseconds

    Returns
    -------
    surf_area : number
    surf_area_err : number
    """
    # TODO
    pass

def clump_dust_luminosity():
    pass

def clump_line_luminosity():
    pass

def compute_physical_conditions(bgps=[], v=2, verbose=True):
    """
    Calculate physical conditions for all clumps. Incudes: dust mass, dust
    luminosity, diameter, dust mass surface density, and surface area.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS catalog to merge data to
    v : number, default 2
        BGPS version, 1 or 2
    verbose : Boolean, default False

    Returns
    -------
    bgps : pd.DataFrame
    """
    if len(bgps) == 0:
        bgps = catalog.read_bgps(exten='all', v=v)
    # new columns
    new_cols = ['dpdf_f', 'dist', 'dist_err', 'mdust', 'mdust_err']
    for col in new_cols:
        bgps[col] = _np.nan
    # dpdf properties
    dpdf = catalog.read_dpdf()
    dpdf_props = _pd.DataFrame(dpdf[1].data)
    dpdf_post = dpdf[5].data
    dist_axis = _np.arange(1000) * 20. + 20.
    lims = [dist_axis.min(), dist_axis.max(), 0, 1]
    for i, cnum in enumerate(dpdf_props.CNUM):
        # retrieve values from catalog
        c_index = _np.argwhere(bgps.cnum == cnum)[0][0]
        c_dpdf = dpdf_post[i]
        flux = bgps.ix[c_index, 'flux']
        # calculate
        dists, dist_probs = mc_sampler_2d(dist_axis, c_dpdf, lims=lims,
            nsample=1e4)
        mdust, mdust_err = clump_simple_dust_mass(dists, flux, tkin=20.)
        # store physical properties
        bgps.ix[c_index, 'dpdf_f'] = 1
        bgps.ix[c_index, 'dist'] = _np.mean(dists)
        bgps.ix[c_index, 'dist_err'] = _np.std(dists)
        bgps.ix[c_index, 'mdust'] = mdust
        bgps.ix[c_index, 'mdust_err'] = mdust_err
        if verbose:
            print '-- clump {} ({:0>3d} / 609)'.format(cnum, i+1)
    return bgps

def plot_dpdf_sampling(n=200):
    """
    Plot a Monte Carlo sampling of a DPDF
    """
    dpdf = catalog.read_dpdf()
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

def print_properties(bgps, out_filen='bgps_props.txt'):
    starless = bgps[(bgps.h2o_f == 0) & (bgps.corn_n == 0) & (bgps.ir_f == 0)]
    h2o_no = bgps[bgps.h2o_f == 0]
    ir_yes = bgps[bgps.ir_f == 1]
    h2o_yes = bgps[bgps.h2o_f == 1]
    hii_yes = bgps[bgps.corn_n > 0]
    ego_yes = bgps[bgps.ego_f > 0]
    df_list = [starless, h2o_no, ir_yes, h2o_yes, hii_yes, ego_yes]
    df_names = ['Starless', 'H2O No', 'IR Yes', 'H2O Yes', 'HII Yes', 'EGO Yes']
    col_list = ['flux', 'h2o_tpk', 'h2o_int', 'h2o_vsp', 'h2o_num_lines']
    for i, df in enumerate(df_list):
        df['nnh/hco'] = df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           ((df.nnh_f == 1) | (df.nnh_f == 3))]['nnh_int'] / \
                        df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           ((df.nnh_f == 1) | (df.nnh_f == 3))]['hco_int']
        df['hco/nh3'] = df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['hco_int'] / \
                        df[((df.hco_f == 1) | (df.hco_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['nh3_pk11']
        df['nnh/nh3'] = df[((df.nnh_f == 1) | (df.nnh_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['hco_int'] / \
                        df[((df.nnh_f == 1) | (df.nnh_f == 3)) &
                           (df.nh3_pk11 / df.nh3_noise11 > 4)]['nh3_pk11']
        # group numbers
        print '-- {}:'.format(df_names[i])
        print 'Number in group: {}'.format(
            df.shape[0])
        print 'Number with DPDFs: {}'.format(
            df[df.dpdf_f > 0].shape[0])
        print 'Number with Tkin: {}'.format(
            df[(df.amm_f > 0) | (df.nh3_mult_n > 0)].shape[0])
        print 'Number with DPDF & Tkin: {}'.format(
            df[((df.amm_f > 0) | (df.nh3_mult_n > 0)) & (df.dpdf_f > 0)].shape[0])
        for col in ['flux', 'flux_40', 'flux_80', 'flux_120']:
            print '{0}: {1} ({2})'.format(col,
                df[col].median(), df[col].shape[0])
        for col in ['hco_tpk', 'hco_int', 'hco_fwhm']:
            print '{0}: {1} ({2})'.format(col,
                df[(df.hco_f == 1) | (df.hco_f == 3)][col].median(),
                df[(df.hco_f == 1) | (df.hco_f == 3)][col].shape[0])
        for col in ['nnh_tpk', 'nnh_int', 'nnh_fwhm']:
            print '{0}: {1} ({2})'.format(col,
                df[(df.nnh_f == 1) | (df.nnh_f == 3)][col].median(),
                df[(df.nnh_f == 1) | (df.nnh_f == 3)][col].shape[0])
        for col in ['nnh/hco', 'hco/nh3', 'nnh/nh3']:
            print '{0}: {1} ({2})'.format(col,
                df[col].median(), df[col].shape[0])
        for col in ['h2o_tpk', 'h2o_int', 'h2o_vsp', 'h2o_num_lines']:
            print '{0}: {1} ({2})'.format(col,
                df[df.h2o_gbt_f > 0][col].median(), df[col].shape[0])
        for col in ['nh3_tkin', 'nh3_pk11']:
            print '{0}: {1} ({2})'.format(col,
                df[df.nh3_pk11 / df.nh3_noise11 > 4][col].median(),
                df[df.nh3_pk11 / df.nh3_noise11 > 4][col].shape[0])
    pass

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
    dpdf = catalog.read_dpdf(v=v)
    flags = _pd.DataFrame(dpdf[1].data)
    for i, row in enumerate(dpdf[5].data):
        _np.savetxt(out_dir + '/v{0}_{1:0>4d}.txt'.format(v,
            flags.CNUM.iloc[i]), row, fmt='%.10e')
    return

def match_dpdf_to_tim_calc(bgps=[]):
    """
    Match BGPS catalog to Tim's calculated quantities.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS dataframe

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read bgps
    if len(bgps) == 0:
        bgps = catalog.read_bgps(exten='all')
    bgps_201 = catalog.read_bgps(v=201)
    bgps_201 = bgps_201[['cnum', 'name']]
    bgps_201 = bgps_201.rename(columns={'cnum': 'v201cnum'})
    # read emaf
    emaf = catalog.read_emaf_dist()
    emaf_cols = ['Seq', 'KDAR', 'dML', 'dMLp', 'dMLm', 'dBAR', 'e_dBAR']
    emaf = emaf[emaf_cols]
    # merge new cnums to bgps
    bgps = _pd.merge(bgps, bgps_201, on='name')
    # merge new emaf to bgps
    bgps = _pd.merge(bgps, emaf, left_on='v201cnum', right_on='Seq', how='outer')
    bgps = bgps.drop(labels=['Seq'], axis=1)
    return bgps
