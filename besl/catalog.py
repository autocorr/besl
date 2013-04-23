"""
Library to match catalogs.
"""
# TODO add survey boundaries to put in null values if not in overlap region
# TODO add in other WISE catalogs
# TODO add merge functions

import os as _os
import sys as _sys
import numpy as _np
import pandas as _pd
import ephem as _ephem
import pywcs as _pywcs
import pyfits as _pyfits

### Directories, paths, and environment variables
class Dirs(object):
    """
    Object to hold directories for interactive editing of paths.
    """
    def __init__(self):
        self.root_dir = '/mnt/eld_data/'
        self.cat_dir = self.root_dir + 'Catalogs/'
        self.bgps_dir = self.root_dir + 'BGPS/Images/v2.0.0/'
        self.working_dir = _os.getcwd() + '/'
        self.out_dir = self.working_dir + 'matched_cats/'
        self.bgps1_filen = 'bgps/bgps_v1.0.csv'
        self.bgps2_filen = 'bgps/bgps_v2.0.csv'
        self.bgps_bounds_filen = 'bgps/bgps_v2.0_bounds.csv'
        self.molcat_filen = 'bgps/bgps_molcat.csv'
        self.wise_filen = 'wise/wise_0-90.csv'
        self.msx_filen = 'red_msx/rms_msx_urquhart.csv'
        self.robit_filen = 'red_robitaille/red_robitaille.csv'
        self.ego_filen = 'ego/ego_{}_cyganowski.csv'
        self.mmb_filen = 'mmb/mmb_all.csv'
        self.gbt_h2o_filen = 'bgps/gbt_h2o_all.csv'
        self.rms_h2o_det_filen = 'red_msx/rms_h2o_det_urquhart.csv'
        self.rms_h2o_noise_filen = 'red_msx/rms_h2o_noise_urquhart.csv'
        self.arcetri_ces_filen = 'arcetri/arcetri_cesaroni.csv'
        self.arcetri_val_filen = 'arcetri/arcetri_valdettaro.csv'
        self.hops_filen = 'hops/hops_walsh.csv'
        self.hrds_ao_filen = 'hrds/hrds_arecibo_slist.csv'
        self.hrds_gbt_filen = 'hrds/hrds_gbt_slist.csv'
        self.hii_bania_filen = 'hrds/known_hii_bania.csv'
        self.v1_dpdf_filen = 'emaf/BGPS_V1_dpdf_table.fits'
        self.v2_dpdf_filen = 'emaf/BGPS_V2_dpdf_table.fits'
        self.v1_emaf_filen = 'emaf/emaf_dist_V1.csv'
        self.v2_emaf_filen = 'emaf/emaf_dist_V2.csv'
        self.gbt_nh3_1_filen = 'bgps/nh3_11B48_fit_objects.sav'
        self.gbt_nh3_2_filen = 'bgps/nh3_bgps_fit_objects.sav'
        self.gbt_nh3_3_filen = 'bgps/nh3_rms_fit_objects.sav'
d = Dirs()

### Read functions ###
# Functions to read catalogs and return pandas DataFrame objects
def read_bgps(v=2):
    """
    Read BGPS catalog, defaults to version 2.0.1. Citation: Ginsburg et al.
    (2013) for v2, Aguirre et al. (2011) for v1.

    Parameters
    ----------
    v : number, 1 or 2
        Catalog version

    Returns
    -------
    bgps : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    if v not in [1, 2]:
        raise ValueError('v = {}. Must be 1 or 2.'.format(v))
    if v == 1:
        bgps = _pd.read_csv(d.cat_dir + d.bgps1_filen, comment='#',
            na_values=['null'], skiprows=4)
        pass
    if v == 2:
        bgps = _pd.read_csv(d.cat_dir + d.bgps2_filen, comment='#',
            na_values=['null'], skiprows=4)
        bgps['cnum'] = _np.arange(1, bgps.shape[0] + 1)
    return bgps

def read_bgps_bounds():
    """
    Read BGPS coordinate boundaries of image files. Citation: Ginsburg et al.
    (2013).

    Returns
    -------
    bgps_bounds : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    bgps_bounds = _pd.read_csv(d.cat_dir + d.bgps_bounds_filen)
    return bgps_bounds

def read_molcat():
    """
    Read BGPS HHT molecular line catalog. Citation: Shirley et al. (2013).

    Returns
    -------
    molcat : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    molcat = _pd.read_csv(d.cat_dir + d.molcat_filen, na_values=['99.99'],
        skiprows=43)
    return molcat

def read_gbt_nh3():
    """
    Read BGPS GBT NH3 observation and temperature fit catalog. Citation:
    Dunham et al. (2011), Rosolowsky et al. (in prep.).

    Returns
    -------
    gbt_nh3 : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    idl_list : list
        IDL save file products
    """
    import idlsave
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nh3_1 = idlsave.read(d.cat_dir + d.gbt_nh3_1_filen, verbose=False)
        nh3_2 = idlsave.read(d.cat_dir + d.gbt_nh3_2_filen, verbose=False)
        nh3_3 = idlsave.read(d.cat_dir + d.gbt_nh3_3_filen, verbose=False)
    idl_list = [nh3_1, nh3_2, nh3_3]
    nh3_1 = _pd.DataFrame(nh3_1.s)
    nh3_2 = _pd.DataFrame(nh3_2.s)
    nh3_3 = _pd.DataFrame(nh3_3.s)
    df_list = [nh3_1, nh3_2, nh3_3]
    gbt_nh3 = _pd.concat(df_list, ignore_index=True)
    gbt_nh3['SNR11'] = gbt_nh3['PK11'] / gbt_nh3['NOISE11']
    gbt_nh3['SNR22'] = gbt_nh3['PK22'] / gbt_nh3['NOISE22']
    gbt_nh3['SNR33'] = gbt_nh3['PK33'] / gbt_nh3['NOISE33']
    return gbt_nh3, idl_list

def read_wise():
    """
    Read WISE catalog. Citation: Wright et al. (2010).

    Returns
    -------
    wise : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    wise = _pd.read_csv(d.cat_dir + d.wise_filen, sep=';', skipinitialspace=True,
        skiprows=116, na_values=['null'])
    # Compute colors
    wise['w1mw2'] = wise.w1mpro - wise.w2mpro
    wise['w1mw3'] = wise.w1mpro - wise.w3mpro
    wise['w1mw4'] = wise.w1mpro - wise.w4mpro
    wise['w2mw3'] = wise.w2mpro - wise.w3mpro
    wise['w2mw4'] = wise.w2mpro - wise.w4mpro
    wise['w3mw4'] = wise.w3mpro - wise.w4mpro
    return wise

def read_wise_sig():
    """
    Read WISE catalog cut for all bands SNR above 7 sigma and no artifact
    flags. Citation: Wright et al. (2010).

    Returns
    -------
    wise : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    wise = read_wise()
    wise = wise[(wise.w1snr > 7) &
                (wise.w2snr > 7) &
                (wise.w3snr > 7) &
                (wise.w4snr > 7)]
    wise = wise[wise.cc_flags == '0000']
    return wise

def read_red_wise(print_diag=False):
    """
    Read WISE catalog and compute colors for YSO candidates. Citation: Wright
    et al. (2010), Koenig et al. (2012).

    Parameters
    ----------
    print_diag : Boolean, default False
        Print diagnostics of catalog when computing

    Returns
    -------
    wise : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    wise = read_wise_sig()
    significant_sources = wise.shape[0]
    # Cut for star-forming galaxies
    # From Koenig et al. (2012) A.1
    wise = wise[(wise.w1mw2 > 0.46 * (wise.w2mw3 - 1.7)) &
                (wise.w1mw2 < -0.06 * (wise.w2mw3 - 4.67)) &
                (wise.w1mw2 > -1.0 * (wise.w2mw3 - 5.1)) &
                (wise.w1mw2 < -4.1 * (wise.w2mw3 - 4.1)) &
                (wise.w2mpro < 12.0) &
                (wise.w2mw3 < 2.3)]
    sfg_sources = significant_sources - wise.shape[0]
    # Cut for broad-line AGN
    # From Koenig et al. (2012) A.2
    wise = wise[((wise.w2mpro < 1.9 * (wise.w2mw3 + 3.16)) &
                 (wise.w2mpro < -1.4 * (wise.w2mw3 - 11.93)) &
                 (wise.w2mpro < 13.5)) |
                ((wise.w1mpro < 1.9 * (wise.w1mw3 + 2.55)) &
                 (wise.w1mpro < 14.0))]
    agn_sources = sfg_sources - wise.shape[0]
    # Cut for shock / PAH emission
    wise = wise[((wise.w1mw2 < 1.0) &
                 (wise.w2mw3 > 2.0)) |
                (((wise.w1mw2 > 1.0) &
                  (wise.w2mw3 < 4.9)) &
                 ((wise.w1mw2 > 0.25) &
                  (wise.w2mw3 < 4.75)))]
    pah_sources = agn_sources - wise.shape[0]
    # From Koenig et al. (2012) A.3
    # TODO Cut for AGB stars
    end_sources = wise.shape[0]
    # Print diagnostic values
    if print_diag is True:
        print '-- Total significant: {}'.format(significant_sources)
        print '-- SFG cut: {}'.format(sfg_sources)
        print '-- AGN cut: {}'.format(agn_sources)
        print '-- PAH cut: {}'.format(pah_sources)
        #print '-- AGB cut: {}'.format(agb_sources)
        print '-- Total after cuts: {}'.format(end_sources)
    return wise

def read_msx():
    """
    Read Red MSX catalog. Citation: Urquhart et al. (2008)

    Returns
    -------
    msx : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    from besl.coord import sexstr2dec as _sexstr2dec
    msx = _pd.read_csv(d.cat_dir + d.msx_filen)
    msx['ra'] = msx['ra_str'].apply(_sexstr2dec, hd='h')
    msx['dec'] = msx['dec_str'].apply(_sexstr2dec, hd='d')
    return msx

def read_robitaille():
    """
    Read Spitzer Red Sources from GLIMPSE. Citation: Robitaille et al. (2008).

    Returns
    -------
    robit : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    robit = _pd.read_csv(d.cat_dir + d.robit_filen, sep=';', skiprows=56,
        skipinitialspace=True)
    return robit

def read_ego():
    """
    Read EGO catalog. Citation: Cyganowski et al. (2008).

    Returns
    -------
    ego : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    df_list : list
        List grouping of individual DataFrames for each
        EGO table (1 - 5)
    """
    skip_list = [91, 15, 30, 15, 19]
    df_list = []
    for i in range(1,6):
        df_list.append(_pd.read_csv(d.cat_dir + d.ego_filen.format(i),
            sep=';', skipinitialspace=True, skiprows=skip_list[i-1]))
    ego = _pd.concat(df_list, ignore_index=True)
    return ego, df_list

def read_mmb():
    """
    Read Methanol Multi-Beam survey catalog. Citation: Green et al. (2010, 2012),
    Caswell et al. (2010, 2011).

    Returns
    -------
    mmb : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    mmb = _pd.read_csv(d.cat_dir + d.mmb_filen, sep=',', na_values=[-999])
    mmb['ra'] = (mmb.ra_h + mmb.ra_m / 60. + mmb.ra_s / (60. * 60.)) \
        * 360. / 24.
    mmb['dec'] = mmb.dec_d + mmb.dec_m / 60. + mmb.dec_s / (60. * 60.)
    return mmb

def read_gbt_h2o():
    """
    Read BGPS GBT H2O maser survey catalog. Citation: Svoboda et al. (2013).

    Returns
    -------
    gbt_h2o : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    tmb2jy = 0.463
    gbt_h2o = _pd.read_csv(d.cat_dir + d.gbt_h2o_filen, na_values=['-999'],
        skiprows=21)
    gbt_h2o['h2o_tpk_jy'] = tmb2jy * gbt_h2o['h2o_tpk']
    gbt_h2o['h2o_tpk_err_jy'] = tmb2jy * gbt_h2o['h2o_tpk_err']
    gbt_h2o['h2o_int_jy'] = tmb2jy * gbt_h2o['h2o_int']
    gbt_h2o['h2o_int_err_jy'] = tmb2jy * gbt_h2o['h2o_int_err']
    return gbt_h2o

def read_rms_h2o():
    """
    Read Red MSX H2O maser catalog. Citation: Urquhart et al. (2009)

    Returns
    -------
    rms_h2o : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    k2jy = 6.41
    rms_h2o = _pd.read_csv(d.cat_dir + d.rms_h2o_det_filen, sep=';',
        skipinitialspace=True, skiprows=47)
    rms_h2o['vsp'] = _np.abs(rms_h2o['Vmin'] - rms_h2o['Vmax'])
    rms_h2o_noise = _pd.read_csv(d.cat_dir + d.rms_h2o_noise_filen,
        sep=';', skipinitialspace=True, skiprows=40)
    rms_h2o_noise['rms_Jy'] = 1e-3 * k2jy * rms_h2o_noise['rms_K']
    rms_h2o = rms_h2o.merge(rms_h2o_noise, how='outer', on='recno')
    rms_h2o['h2o_f'] = 1
    rms_h2o['h2o_f'][rms_h2o.Fpeak.isnull()] = 0
    return rms_h2o

def read_arcetri_cesaroni():
    """
    Read Arcetri H2O maser survey catalog. Citation: Cesaroni et al. 1988.
    Note, positions are in 1950 epoch coordinates.

    Returns
    -------
    arcetri_ces : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    from besl.coord import sexstr2dec as _sexstr2dec
    arcetri_ces = _pd.read_csv(d.cat_dir + d.arcetri_ces_filen, sep=';',
        skipinitialspace=True, skiprows=1)
    arcetri_ces['ra'] = arcetri_ces['ra_str'].apply(_sexstr2dec, hd='h')
    arcetri_ces['dec'] = arcetri_ces['dec_str'].apply(_sexstr2dec, hd='d')
    return arcetri_ces

def read_arcetri_valdettaro():
    """
    Read Arcetri H2O maser survey catalog. Citation: Valdettaro et al. (2004).

    Returns
    -------
    arcetri_val : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    arcetri_val = _pd.read_csv(d.cat_dir + d.arcetri_val_filen, sep=';',
        skipinitialspace=True, skiprows=71)
    arcetri_val['vsp'] = _np.abs(arcetri_val['Vmin'] - arcetri_val['Vmax'])
    arcetri_val['snr'] = arcetri_val['Fp'] / arcetri_val['Sig']
    arcetri_val['h2o_f'] = 0
    arcetri_val['h2o_f'][arcetri_val.snr > 3] = 1
    return arcetri_val

def read_hops():
    """
    Read H2O Southern Galactic Plane Survey (HOPS) catalog. Citation: Walsh et
    al. (2012).

    Returns
    -------
    hops : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    k2jy = 12.5
    hops = _pd.read_csv(d.cat_dir + d.hops_filen)
    hops['Tbdv_Jykms'] = k2jy * hops['Tbdv_Kkms']
    hops['dTbdv_Jykms'] = k2jy * hops['dTbdv_Kkms']
    hops['T_peak_Jy'] = k2jy * hops['T_peak_K']
    hops['RMS_Jy'] = k2jy * hops['RMS_K']
    hops['vsp'] = _np.abs(hops['vMin_kms'] - hops['vMax_kms'])
    return hops

def read_hrds_arecibo():
    """
    Read HII Region Discovery Survey (HRDS) Arcibo catalog. Citation: Bania et
    al. ().

    Returns
    -------
    hrds_ao : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    hrds_ao = _pd.read_csv(d.cat_dir + d.hrds_ao_filen, sep=';',
        skipinitialspace=True, skiprows=36)
    return hrds_ao

def read_hrds_gbt():
    """
    Read HII Region Discovery Survey (HRDS) GBT catalog. Citation: Anderson et
    al. ().

    Returns
    -------
    hrds_gbt : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    hrds_gbt = _pd.read_csv(d.cat_dir + d.hrds_gbt_filen, sep=';',
        skipinitialspace=True, skiprows=49)
    return hrds_gbt

def read_hii_bania():
    """
    Read known HII region catalog. Citation: Bania et al. ().

    Returns
    -------
    hii_bania : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    hii_bania = _pd.read_csv(d.cat_dir + d.hii_bania_filen, sep=';',
        skipinitialspace=True, skiprows=8)
    return hii_bania

def read_dpdf(v=2):
    """
    Read Distance Probability Distribution Functions. Citation:
    Ellsworth-Bowers et al. (2013).

    Parameters
    ----------
    v : number {1, 2}, default 2
        Version of BGPS to use

    Returns
    -------
    dpdf : pandas.DataFrame
    """
    if v not in [1, 2]:
        raise ValueError('v = {}. Must be 1 or 2.'.format(v))
    if v == 1:
        dpdf = _pyfits.open(d.cat_dir + d.v1_dpdf_filen)
        return dpdf
    if v == 2:
        dpdf = _pyfits.open(d.cat_dir + d.v2_dpdf_filen)
        return dpdf

def read_emaf_dist(v=2):
    """
    Read EMAF calculated distances and parameters. Citation: Ellsworth-Bowers
    et al. (2013).

    Parameters
    ----------
    v : number {1, 2}, default 2
        Version of BGPS to use

    Returns
    -------
    emaf : pandas.DataFrame
    """
    if v not in [1, 2]:
        raise ValueError('v = {}. Must be 1 or 2.'.format(v))
    if v == 1:
        emaf = _pd.read_csv(d.cat_dir + d.v1_emaf_filen, sep=';',
            skipinitialspace=True, skiprows=33, na_values=['---'])
        return emaf
    if v == 2:
        emaf = _pd.read_csv(d.cat_dir + d.v2_emaf_filen, sep=';',
            skipinitialspace=True, skiprows=31, na_values=['---'])
        return emaf

def select_bgps_field(lon, lat, coord_type='eq'):
    """
    Select BGPS image field name from sky position in galactic longitude and
    latitude.

    Parameters
    ----------
    glon : number
        Galactic longitude in decimal degrees
    glat : number
        Galactic latitude in decimal degrees

    Returns
    -------
    field : string
        Image field name
    """
    from besl.coord import eq2gal
    if coord_type not in ['eq', 'gal']:
        raise ValueError(
            'coord_type = {}. Must be eq or gal.'.format(coord_type))
    # read catalog
    bgps_bounds = read_bgps_bounds()
    # fudge factor because of frame spacing
    epsilon = 0.0015
    bgps_bounds['glon_min'] -= epsilon
    bgps_bounds['glon_max'] += epsilon
    # convert equatorial to galactic
    if coord_type == 'eq':
        lon, lat = eq2gal(lon, lat)
    # select field
    glon_max_sep = _np.mod(lon - bgps_bounds.glon_max + 180, 360) - 180
    glon_min_sep = _np.mod(lon - bgps_bounds.glon_min + 180, 360) - 180
    field = bgps_bounds.ix[_np.argwhere((glon_max_sep.values > 0) &
        (glon_min_sep.values < 0))[0][0]]['field']
    #field = bgps_bounds[(bgps_bounds.glon_min < lon) &
    #                    (bgps_bounds.glon_max > lon) &
    #                    (bgps_bounds.glat_min < lat) &
    #                    (bgps_bounds.glat_max > lat)].field
    if len(field) == 0:
        return _np.nan
    else:
        return field

### Clump matching
# Procedures dealing with the BGPS label masks
def clump_match(haystack_list, cnum, coord_type='eq', pix_size=7.5, bgps_ver=2):
    """
    Match a catalog of coordinates to the coordinates found in the pixels of
    the label mask for a BGPS clump. Returns a list of indices of catalog
    sources that are within the label mask. Match radius circumscribes square
    pixel, ie a circle centered on the pixel center with radius drawn from the
    center to a corner (sqrt(2) * pix_size / 2).

    Parameters
    ----------
    haystack_list : array like
        List of coordinates in (RA, Dec) or (Lon, Lat) along axis = 0. Values
        must be in decimal degrees.
    cnum : number
        Catalog number of BGPS clump
    coord_type : string ('eq', 'gal'), default 'eq'
        Celestial coordinate type, select either 'eq' = equatorial or 'gal' =
        Galactic.
    pix_size : number, default 7.5"
        Pixel size in arcseconds
    bgps_ver : number (1, 2), default 2
        Version of BGPS, v1.0.1 or v2.0.1.

    Returns
    -------
    match_index : list
        List of indices for sources that match within label mask
    """
    # TODO support v1
    import ipdb as pdb
    if bgps_ver not in [1, 2]:
        raise ValueError('bgps_ver = {}. Must be 1 or 2.'.format(v))
    from besl.coord import nearest_match_coords
    search_rad = _np.sqrt(2) * pix_size / 2.
    index_list = []
    # read in catalog
    bgps_bounds = read_bgps_bounds()
    bgps = read_bgps(v=bgps_ver)
    field = bgps[bgps.cnum == cnum]['field'].values[0]
    fcnum = bgps[bgps.cnum == cnum]['field_cnum'].values[0]
    # read in label mask
    bgps_label = _pyfits.open(d.bgps_dir +
        'v2.0_ds2_{}_13pca_labelmask.fits'.format(field))
    bgps_wcs = _pywcs.WCS(bgps_label[0].header)
    # select coordinates of all pixels that match
    matched_pixels = _np.argwhere(bgps_label[0].data == fcnum)
    pixel_coords = bgps_wcs.all_pix2sky(matched_pixels[:,1],
        matched_pixels[:,0], 0) # origin = 0 for numpy
    pixel_coords = _np.transpose(pixel_coords)
    # for each pixel, use as "needle" in the haystack of the supplied catalog
    for coord in pixel_coords:
        # if haystack in Equatorial, convert pixel from Galactic
        if coord_type == 'eq':
            gal = _ephem.Galactic(_np.deg2rad(coord[0]), _np.deg2rad(coord[1]),
                epoch='2000')
            coord = _np.degrees(gal.to_radec())
        (matchn, min_index, min_dist) = nearest_match_coords(coord,
            haystack_list, min_sep=search_rad)
        if matchn != 0:
            index_list.append(min_index)
    index_list = _np.unique(index_list)
    return index_list

### Merge functions to match clumps
def match_h2o_maser():
    # TODO
    pass

def match_mmb():
    # TODO
    pass

def match_wise():
    # TODO
    pass

def match_IR():
    # TODO
    pass

def match_all_clumps():
    # TODO
    # select clumps in overlap range or put null value if not in overlap range
    # loop through relevant clumps
    #   select only clump in label mask fits
    pass

### DS9 regions
def create_point_region(lon, lat, out_filen='ds9.reg', marker_type='circle',
        coord_type='fk5', color='green'):
    """
    Create a DS9 region file from a list of longitude and latitude coordinates.

    Parameters
    ----------
    lon : array
    lat : array
    out_filen : string, default 'ds9.reg'
        Filename of output DS9 region file.
    marker_type : string, default 'circle'
        Region marker type. Supported: {circle, box, diamond, cross, x, arrow,
        boxcircle}.
    coord_type : string, default 'fk5'
        Celestial coordinate type. Supported: {image, linear, fk4, fk5,
        galactic, ecliptic, icrs, physical, amplifier, detector}.
    color : string, default 'green'
        Region color. Supported: {white, black, red, green, blue, cyan, magenta,
        yellow}.

    Returns
    -------
    out_file : regions file
    """
    point_strings = ['circle', 'box', 'diamond', 'cross', 'x', 'arrow',
        'boxcircle']
    out_file = open(d.out_dir + out_filen, 'w')
    out_file.write('# global color={0}'.format(color))
    out_file.write('{0}'.format(coord_type))
    if marker in point_types:
        for i in xrange(lon.shape[0]):
            out_file.write('{marker} point {lon} {lat}'.format(
                marker=marker,
                lon=lon[i],
                lat=lat[i]))
    out_file.close()
    print '-- DS9 region file written to {}'.format(out_filen)
    return

