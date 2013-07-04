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
import catalog, units
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

def mc_sampler_1d(x, y, lims=[0,1,0,1], nsample=1e3):
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

def tkin_distrib(tkins):
    """
    Interpolate the distribution of a series of temperatures.

    Parameters
    ----------
    tkins : list-like
        Series of temperatures

    Returns
    -------
    tkin_fn : scipy.interpolate.interpolate.interp1d
    """
    tkins = _np.sort(tkins)
    kernel = gaussian_kde(tkins)
    tkin_fn = interp1d(tkins, kernel(tkins))
    return tkin_fn

def draw_dpdf_samples(cnum, nsample=1e4):
    """
    Draw a set of sampled distances for a BGPS clump with a DPDF.

    Parameters
    ----------
    cnum : number
        BGPS v2 catalog number

    Returns
    -------
    dist_draws : np.array
    """
    # Read in DPDFs
    dpdf = catalog.read_dpdf()
    # Determine catalog index number
    flags = _pd.DataFrame(dpdf[1].data)
    if len(flags) == 0:
        raise Exception('cnum not in DPDFs: {0}'.format(cnum))
    i = flags.index[flags.CNUM == cnum][0]
    # Select data
    x = _np.arange(1000) * 20. + 20.
    y = dpdf[5].data[i]
    lims = [x.min(), x.max(), 0, 1]
    # Draw distance samples
    dist_draws, dist_probs = mc_sampler_1d(x, y, lims=lims, nsample=nsample)
    return dist_draws

def generate_mass_samples(cnum, tkin_err=None):
    """
    Generate a sample of masses with values drawn from the DPDFs and kinetic
    temperature distributions.
    """
    # TODO
    # draw dml samples
    # draw tkin samples
    #   use tkin error and normal deviates if given value
    if tkin_err is None:
        tkin_draw = np.random.normal(loc=tkin, scale=tkin_err, size=1e4)
    # normalize output distributio
    # return both a kde and hist
    return

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

def clump_simple_dust_mass(dist, snu11, tkin=20., use_kpc=False):
    """
    Calculate the dust mass in Solar masses based on an OH5 dust opacity model
    at 1.1mm, as used by the BGPS.

    Parameters
    ----------
    dist : number
        Distance in pc
    snu11 : number
        1.1 mm flux density in Jy
    tkin : number, default 20.
        Kinetic Temperature in K
    use_kpc : Bool, default False
        Use kpc instead of pc

    Returns
    -------
    mdust : number
    """
    if use_kpc:
        dist_factor = 1
    else:
        dist_factor = 1e-3
    mdust = 14.26 * (_np.exp(13.01 / tkin) - 1) * snu11 * (dist *
        dist_factor)**2
    return mdust

def clump_surface_area(dist, area, use_kpc=False):
    """
    Calculate the clump surface area in square pc based on the BGPS label mask
    area.

    Parameters
    ----------
    dist : number
        distance in pc
    area : number
        Surface area in square arcsec
    use_kpc : Bool, default False
        Use kpc instead of pc

    Returns
    -------
    area_dist : number
        Surface area at distance in square pc
    """
    if use_kpc:
        dist_factor = 1
    else:
        dist_factor = 1e3
    area_dist = (units.cgs.au / units.cgs.pc)**2 * (dist * dist_factor)**2 * \
        area
    return area_dist

def clump_diameter(dist, area, use_kpc=False):
    """
    Calculate the clump diameter in pc based on the BGPS label mask area.

    Parameters
    ----------
    dist : number
        distance in pc
    area : number
        Surface area in square arcsec
    use_kpc : boolean, default False
        Use kpc instead of pc

    Returns
    -------
    diam_dist : number
        diameter at distance in pc
    """
    if use_kpc:
        dist_factor = 1
    else:
        dist_factor = 1e3
    diam_dist = (units.cgs.au / units.cgs.pc) * (dist * dist_factor) * \
        2 * _np.sqrt(area / _np.pi)
    return diam_dist

def calc_ml_physical_conditions(bgps=[], neighbor=False):
    """
    Calculate the physical conditions for all clumps using the maximum
    likelihood distance. Includes: dust mass, diameter, and surface area.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS dataframe, if left blank, will read exten='all'
    neighbor : boolean
        Include neighbor KDAR resolutions in calculations

    Returns
    -------
    bgps : pd.DataFrame
    """
    if len(bgps) == 0:
        bgps = catalog.read_bgps(exten='all')
    if 'dpdf_dML' not in bgps.columns:
        raise ValueError('Incorrect columns')
    good_KDARs = ['N', 'F', 'T']
    if neighbor:
        bgps['all_dML'] = bgps['neighbor_dML'].combine_first(bgps['dpdf_dML'])
        dML_col = 'all_dML'
    else:
        dML_col = 'dpdf_dML'
    bgps['dust_mass'] = bgps[[dML_col, 'flux', 'nh3_tkin']].apply(lambda row:
        clump_simple_dust_mass(row[dML_col], row['flux'], row['nh3_tkin'],
        use_kpc=True), axis=1)
    bgps['avg_diam'] = bgps[[dML_col, 'rind_area']].apply(lambda row:
        clump_diameter(row[dML_col], row['rind_area'], use_kpc=True),
        axis=1)
    bgps['rind_surf_area'] = bgps[[dML_col, 'rind_area']].apply(lambda row:
        clump_surface_area(row[dML_col], row['rind_area'], use_kpc=True),
        axis=1)
    return bgps

def calc_physical_conditions(bgps=[], v=2, verbose=True):
    """
    Calculate physical conditions for all clumps. Incudes: dust mass, dust
    luminosity, diameter, dust mass surface density, and surface area.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS catalog to merge data to
    v : number, default 2
        BGPS version, 1 or 2
    verbose : boolean, default False

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
        dists, dist_probs = mc_sampler_1d(dist_axis, c_dpdf, lims=lims,
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

def gen_stages(bgps=[], stages_group=2, label=None):
    """
    Generate a list of stages from the BGPS catalog.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS catalog in a pandas dataframe, requires evolutionary flags from an
        all-matched catalog.
    stages_group : number, default 3
        Divisions for grouped stages in return list. Options include:
            * 0 (2) : starless, protostellar
            * 1 (3) : starless, IR Y / H2O N, IR Y / H2O Y
            * 2 (6) : starless, H2O N, IR Y, H2O Y, EGO Y, UCHII Y
    label : string
        Replace column with NaN's for log plots

    Returns
    -------
    stages : list
    """
    # Nan values for log-plots
    if label is not None:
        bgps[label][bgps[label] <= 0] = _np.nan
    # evo stages
    if stages_group == 0:
        starless = bgps[(bgps.h2o_f == 0) & (bgps.corn_n == 0) & (bgps.ir_f == 0)]
        protostellar = bgps[(bgps.h2o_f == 1) | (bgps.corn_n > 0) | (bgps.ir_f
            == 1)]
        stages = [starless, protostellar]
        return stages
    elif stages_group == 1:
        starless = bgps[(bgps.h2o_f == 0) & (bgps.corn_n == 0) & (bgps.ir_f == 0)]
        h2o_no = bgps[(bgps.h2o_f == 0 ) & (bgps.ir_f == 1)]
        h2o_yes = bgps[(bgps.h2o_f > 0) & (bgps.ir_f == 1)]
        stages = [starless, h2o_no, h2o_yes]
        return stages
    elif stages_group == 2:
        starless = bgps[(bgps.h2o_f == 0) & (bgps.corn_n == 0) & (bgps.ir_f == 0)]
        h2o_no = bgps[bgps.h2o_f == 0]
        ir_yes = bgps[bgps.ir_f == 1]
        h2o_yes = bgps[bgps.h2o_f == 1]
        ego_yes = bgps[bgps.ego_n > 0]
        hii_yes = bgps[bgps.corn_n > 0]
        stages = [starless, h2o_no, ir_yes, h2o_yes, ego_yes, hii_yes]
        return stages
    else:
        raise ValueError('Invalid stages_group: {0}.'.format(stages_group))

def write_dpdf_outfiles(out_dir='dpdf_ascii', v=2):
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

def write_emaf_table():
    """
    Create a CSV table from DPDF fits table in Ellsworth-Bowers et al. (2013).
    """
    dpdf = catalog.read_dpdf()
    tab1 = _pd.DataFrame(dpdf[1].data)
    tab1 = tab1.rename(columns={key: 'dpdf_' + key.lower() + '_f' for key in
        tab1.columns})
    tab1 = tab1.rename(columns={'dpdf_cnum_f': 'v201cnum'})
    tab2 = []
    tab2_header = ['dpdf_' + elem for elem in ['dML', 'dMLm', 'dMLp', 'dBAR',
    'dBAR_err', 'PML', 'FW68', 'dtan', 'KDAR']]
    for row in dpdf[6].data:
        tab2.append([row[0][0], row[0][1], row[0][2], row[1][0], row[1][1],
            row[3], row[4], row[5], row[6]])
    tab2 = _pd.DataFrame(tab2, columns=tab2_header)
    emaf = tab1.join(tab2)
    emaf.to_csv('emaf_dist_V2.csv', index=False)
    return emaf


