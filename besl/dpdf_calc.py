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
import catalog
import units
from .catalog import read_cat
from .util import convert_endian
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
    if len(tkins) < 20:
        raise ValueError('tkins : {0} < 20.'.format(len(tkins)))
    tkins = _np.sort(tkins)
    kernel = gaussian_kde(tkins)
    tkin_fn = interp1d(tkins, kernel(tkins))
    return tkin_fn


def draw_tkin_samples(tkin_fn=None, tkin=None, tkin_err=None, nsample=1e2):
    """
    Draw a set of sampled NH3 kinetic temperatures from a distribution or use
    Gaussian random deviates if measurement known.

    Parameters
    ----------
    tkin_fn : scipy.interpolate.interpolate.interp1d
        Interpolated function for distribution of temperatures
    tkin : number
        T_K in Kelvin
    tkin_err : number
        Uncertainty in T_K in Kelvin
    nsample : number, default 1e4
        Number of samples to draw

    Returns
    -------
    tkin_draw : np.array
    """
    if tkin_fn is not None:
        x = tkin_fn.x
        y = tkin_fn.y
        lims = [x.min(), x.max(), y.min(), y.max()]
        tkin_draw, tkin_probs = mc_sampler_1d(x, y, lims=lims, nsample=nsample)
        return tkin_draw
    elif (tkin is not None) & (tkin_err is not None):
        tkin_draw = _np.random.normal(loc=tkin, scale=tkin_err, size=nsample)
        return tkin_draw
    else:
        raise Exception('Unexpected exception.')


def draw_dist_samples(cnum, nsample=1e2):
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
    i = flags.index[flags.CNUM == cnum][0]
    # Select data
    x = _np.arange(1000) * 20. + 20.
    y = dpdf[5].data[i]
    lims = [x.min(), x.max(), 0, 1]
    # Draw distance samples
    dist_draws, dist_probs = mc_sampler_1d(x, y, lims=lims, nsample=nsample)
    return dist_draws


def draw_mass_samples(cnum, snu11, tkin_fn=None, tkin=None, tkin_err=None,
    nsample=1e2):
    """
    Generate a sample of masses with values drawn from the DPDFs and kinetic
    temperature distributions.
    """
    # Draw dml samples
    dist_samples = draw_dist_samples(cnum, nsample=nsample)
    # Draw tkin samples
    tkin_samples = draw_tkin_samples(tkin_fn=tkin_fn, tkin=tkin,
        tkin_err=tkin_err, nsample=nsample)
    # Calculate masses
    mass_samples = clump_simple_dust_mass(dist_samples, snu11,
        tkin=tkin_samples)
    return mass_samples


def gen_stage_mass_samples(stage, nsample=1e2):
    """
    Generate a sample of masses for a set of clumps in an evolutionary stage.

    Parameters
    ----------
    stage : pd.DataFrame

    Returns
    -------
    stage_samples : np.array
    """
    good_kdars = ['T', 'N', 'F']
    # Make temperature distribution
    tkins = stage['nh3_tkin'].values
    tkins = tkins[_np.isfinite(tkins)]
    tkin_fn = tkin_distrib(tkins)
    # Monte Carlo sample each clump in stage
    stage_samples = []
    for i in stage.index:
        # Check for DPDF
        kdar = stage.ix[i, 'dpdf_KDAR']
        nkdar = stage.ix[i, 'neighbor_KDAR']
        if kdar in good_kdars:
            cnum = stage.ix[i, 'v201cnum']
        elif nkdar in good_kdars:
            cnum = stage.ix[i, 'neighbor_KDAR_cnum']
        else:
            continue
        # Properties for the clump
        snu11 = stage.ix[i, 'flux']
        tkin = stage.ix[i, 'nh3_tkin']
        tkin_err = stage.ix[i, 'nh3_tkin_err']
        # Use temperature distribution
        if (_np.isnan(tkin)) | (tkin < 7) | (tkin > 150):
            draw = draw_mass_samples(cnum=cnum, snu11=snu11, tkin_fn=tkin_fn,
                nsample=nsample)
        # Use measured temperature and normal deviates
        elif (_np.isfinite(tkin)) & (tkin > 7) & (tkin < 150):
            draw = draw_mass_samples(cnum=cnum, snu11=snu11, tkin=tkin,
                tkin_err=tkin_err, nsample=nsample)
        else:
            raise Exception('Unexpected error.')
        stage_samples.append(draw)
    stage_samples = _np.ravel(stage_samples)
    return stage_samples


def gen_stage_area_samples(stage, nsample=1e2, radius=False, flatten=True):
    """
    Generate a set of area samples from the DPDFs.

    Parameters
    ----------
    stage : pd.DataFrame
        BGPS catalog for a selected evolutionary stage
    nsample : number
        Number of samples to draw for each good DPDF clump
    radius : boolean
        Return mean radius instead of area
    flatten : boolean
        Flatten the output data from a row for each clump's samples or into a
        flat array

    Returns
    -------
    stage_samples : np.array
    """
    good_kdars = ['T', 'N', 'F']
    # Monte Carlo sample each clump in stage
    stage_samples = []
    for i in stage.index:
        area = stage.ix[i, 'rind_area']
        kdar = stage.ix[i, 'dpdf_KDAR']
        nkdar = stage.ix[i, 'neighbor_KDAR']
        # Check for DPDF
        if kdar in good_kdars:
            cnum = stage.ix[i, 'v201cnum']
        elif nkdar in good_kdars:
            cnum = stage.ix[i, 'neighbor_KDAR_cnum']
        else:
            continue
        draw = draw_dist_samples(cnum=cnum, nsample=nsample)
        draw = clump_surface_area(draw, area)
        stage_samples.append(draw)
    stage_samples = _np.array(stage_samples)
    if flatten:
        stage_samples = _np.ravel(stage_samples)
    if radius:
        return _np.sqrt(stage_samples / _np.pi)
    else:
        return stage_samples


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
    bnu = planck_fn(nu, tkin, freq=True)
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


def calc_ml_physical_conditions(bgps=None, neighbor=False):
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
    if bgps is None:
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


def evo_stages(bgps=None, stages_group=2, label=None):
    """
    Generate a list of stages from the BGPS catalog.

    Parameters
    ----------
    bgps : pd.DataFrame
        BGPS catalog in a pandas dataframe, requires evolutionary flags from an
        all-matched catalog.
    stages_group : number, default 2
        Divisions for grouped stages in return list. Options include:
            * 0 (2) : Starless, Protostellar
            * 1 (3) : Starless, IR Y / H2O N, IR Y / H2O Y
            * 2 (7) : Starless, H2O N, IR Y, H2O Y, CH3OH Y, UCHII Y
            * 3 (X) : All survey flags and group flags
    label : string
        Replace negative column values associated to a label string with NaN's
        for log plots.

    Returns
    -------
    stages : pd.DataFrame
        List of stages `pd.DataFrames`
    anno_labels : dict
        Dictionary of label names
    """
    # Nan values for log-plots
    if bgps is None:
        bgps = read_cat('bgps_v210_evo')
    if label is not None:
        bgps[label][bgps[label] <= 0] = _np.nan
    # Evolutionary stage sets
    if stages_group == 0:
        starless = bgps[bgps.sf_f == 0]
        protostellar = bgps[bgps.sf_f == 1]
        stages = [starless, protostellar]
        anno_labels = [
            {'label': 'Starless'},
            {'label': 'Protostellar'}]
        return stages, anno_labels
    elif stages_group == 1:
        starless = bgps[bgps.sf_f == 0]
        h2o_no = bgps[(bgps.h2o_f == 0) & (bgps.ir_f == 1)]
        h2o_yes = bgps[(bgps.h2o_f > 0) & (bgps.ir_f == 1)]
        stages = [starless, h2o_no, h2o_yes]
        anno_labels = [
            {'label': r'${\rm Starless}$'},
            {'label': r'${\rm H_2O \ \  N}$'},
            {'label': r'${\rm H_2O \ \ Y}$'}]
        return stages, anno_labels
    elif stages_group == 2:
        starless = bgps[bgps.sf_f == 0]
        h2o_no = bgps[bgps.h2o_f == 0]
        ir_yes = bgps[bgps.ir_f == 1]
        h2o_yes = bgps[bgps.h2o_f == 1]
        ch3oh_yes = bgps[bgps.ch3oh_f == 1]
        hii_yes = bgps[bgps.corn_n > 0]
        stages = [starless, h2o_no, ir_yes, h2o_yes, ch3oh_yes, hii_yes]
        anno_labels = [
            {'label': r'${\rm Starless}$'},
            {'label': r'${\rm H_2O \ \  N}$'},
            {'label': r'${\rm IR \ \ Y}$'},
            {'label': r'${\rm H_2O \ \ Y}$'},
            {'label': r'${\rm CH_3OH \ \ Y}$'},
            {'label': r'${\rm UCH\sc{II} \ \ Y}$'}]
        return stages, anno_labels
    elif stages_group == 3:
        starless = bgps[bgps.sf_f == 0]
        h2o_no = bgps[bgps.h2o_f == 0]
        ir_yes = bgps[bgps.ir_f == 1]
        robit_yes = bgps[(bgps.robit_n > 0) & (bgps.robit_f > 0)]
        ego_yes = bgps[bgps.ego_n > 0]
        rms_yes = bgps[(bgps.red_msx_n > 0) & (bgps.red_msx_f > 0)]
        h2o_yes = bgps[bgps.h2o_f == 1]
        gbt_yes = bgps[(bgps.h2o_gbt_n > 0) & (bgps.h2o_gbt_f > 0)]
        arc_yes = bgps[(bgps.h2o_arc_n > 0) & (bgps.h2o_arc_f > 0)]
        hops_yes = bgps[(bgps.h2o_hops_n > 0)]
        ch3oh_yes = bgps[bgps.ch3oh_f == 1]
        pandi_yes = bgps[(bgps.ch3oh_pandi_n > 0)]
        mmb_yes = bgps[(bgps.ch3oh_mmb_n > 0)]
        pesta_yes = bgps[(bgps.ch3oh_pesta_n > 0)]
        hii_yes = bgps[bgps.corn_n > 0]
        stages = [starless, h2o_no, ir_yes, robit_yes, ego_yes, rms_yes,
                  h2o_yes, gbt_yes, arc_yes, hops_yes, ch3oh_yes, pandi_yes,
                  mmb_yes, pesta_yes, hii_yes]
        anno_labels = [
            {'label': r'${\rm Starless}$'},
            {'label': r'${\rm H_2O \ \  N}$'},
            {'label': r'${\rm IR \ \ Y}$'},
            {'label': r'${\rm Robitaille \ \ Y}$'},
            {'label': r'${\rm EGO \ \ Y}$'},
            {'label': r'${\rm RMS \ \ Y}$'},
            {'label': r'${\rm H_2O \ \ Y}$'},
            {'label': r'${\rm GBT \ \ Y}$'},
            {'label': r'${\rm Arcetri \ \ Y}$'},
            {'label': r'${\rm HOPS \ \ Y}$'},
            {'label': r'${\rm CH_3OH \ \ Y}$'},
            {'label': r'${\rm MMB \ \ Y}$'},
            {'label': r'${\rm Pandian \ \ Y}$'},
            {'label': r'${\rm Pestalozzi \ \ Y}$'},
            {'label': r'${\rm UCH\sc{II} \ \ Y}$'}]
        return stages, anno_labels
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
    if not _os.path.exists(out_dir):
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
    Only works on v2.1.0.
    """
    # Read in data
    dpdf = catalog.read_dpdf(v=2)
    tab1 = _pd.DataFrame(dpdf[1].data)
    tab1 = tab1.rename(columns={key: 'dpdf_' + key.lower() + '_f' for key in
        tab1.columns})
    tab1 = tab1.rename(columns={'dpdf_cnum_f': 'v210cnum'})
    tab2 = []
    tab2_header = ['dpdf_' + elem for elem in ['dML', 'dMLm', 'dMLp', 'dBAR',
                   'dBAR_err', 'PML', 'FW68', 'dtan', 'KDAR']]
    for row in dpdf[9].data:
        tab2.append([row[0][0], row[0][1], row[0][2], row[1][0], row[1][1],
            row[3], row[4], row[5], row[6]])
    tab2 = _pd.DataFrame(tab2, columns=tab2_header)
    # Table for velocity information
    tab3 = _pd.DataFrame(dpdf[10].data)
    tab3 = tab3.rename(columns={key: 'dpdf_' + key.lower() for key in tab3.columns})
    # Convert to system Endian format
    tab1 = convert_endian(tab1)
    tab2 = convert_endian(tab2)
    tab3 = convert_endian(tab3)
    # Merge tables
    emaf = tab1.join(tab2)
    emaf = emaf.join(tab3)
    emaf = emaf.rename(columns={'dpdf_glon_f': 'dpdf_glon',
                                'dpdf_glat_f': 'dpdf_glat'})
    emaf.to_csv('bgps_v210_dpdf_props.csv', index=False)
    return emaf


