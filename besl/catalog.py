"""
=======================
Read and Match Catalogs
=======================

Library to match catalogs.

"""

import os as _os
import numpy as _np
import pandas as _pd
import ephem as _ephem
from tempfile import TemporaryFile
from astropy import wcs
from astropy.io import (ascii, fits)
from scipy.interpolate import interp1d
from .coord import eq2gal, pd_eq2gal
from .mathf import ang_diff, bin_factor
from .paths import all_paths as d


### Read functions ###
# Functions to read catalogs and return pandas DataFrame objects
def read_cat(filen=None, print_cats=False):
    """
    Read collected catalogs from the `_collected` directory. CSV formats should
    have clean syntax for `pd.read_csv`.

    Parameters
    ----------
    filen : str
        Catalog file name, without extension
    print_cats : bool, default False
        Print available catalogs

    Returns
    -------
    df : pd.DataFrame
        DataFrame of catalog.
    """
    base_path = d.cat_dir + d.collected
    if print_cats:
        print '-- Available catalogs:'
        for f in _os.listdir(base_path):
            f, _ = _os.path.splitext(f)
            print '   {0}'.format(f)
        return
    if filen is None:
        raise ValueError('Catalog name is None.')
    df = _pd.read_csv(base_path + filen + '.csv')
    return df


def read_bgps(exten='none', v=200):
    """
    Read BGPS catalog, defaults to version 2.0.1. Citation: Ginsburg et al.
    (2013) for v2, Aguirre et al. (2011) for v1.

    Parameters
    ----------
    exten : str, default 'none'
        BGPS extension.
        'none' -> default BGPS
        'all'  -> super matched BGPS
        'cols' -> print column descriptions
    v : number, {1, 2, 201}, default 2
        Catalog version

    Returns
    -------
    bgps : pandas.DataFrame
        Output catalog
    """
    vers = {100: '1.0.0',
            101: '1.0.1',
            200: '2.0',
            201: '2.0.1',
            '2d': '2.0d'}
    if v not in vers.keys():
        raise ValueError('Invalid version, v = {}.'.format(v))
    if exten == 'none':
        if (v == 101) | (v == 100):
            bgps = _pd.read_csv(d.cat_dir + d.bgps_filen.format(vers[v]),
                comment='#', na_values=['null'], skiprows=4)
            return bgps
        elif v == 200:
            bgps = _pd.read_csv(d.cat_dir + d.bgps_filen.format(vers[v]),
                comment='#', na_values=['null'], skiprows=4)
            bgps['cnum'] = _np.arange(1, bgps.shape[0] + 1)
            return bgps
        elif v == 201:
            bgps = _pd.read_csv(d.cat_dir + d.bgps_filen.format(vers[v]),
                na_values=['---'])
            return bgps
        elif v == '2d':
            bgps = _pd.read_csv(d.cat_dir + d.bgps_filen.format(vers[v]))
            return bgps
    elif exten == 'all':
        bgps = _pd.read_csv(d.cat_dir + d.bgps_ext_filen.format(vers[2], 'all',
            'csv'))
        return bgps
    elif exten == 'cols':
        cols = open(d.cat_dir + d.bgps2_ext_filen.format(vers[2], 'columns',
            'txt'))
        print cols.read()
        return
    else:
        raise ValueError('Invalid file extension: exten = {}.'.format(exten))


def read_bgps_bounds():
    """
    Read BGPS coordinate boundaries of image files. Citation: Ginsburg et al.
    (2013).

    Returns
    -------
    bgps_bounds : pandas.DataFrame
        Output catalog
    """
    bgps_bounds = _pd.read_csv(d.cat_dir + d.bgps_bounds_filen)
    return bgps_bounds


def read_bgps_vel(no_clobber=False):
    """
    Read the BGPS collected velocity catalog. Contains HCO+, N2H+, CS, NH3, and
    On-Off GRS 13CO. Citation: Ellsworth-Bowers (personal communication).

    Parameters
    ----------
    no_clobber : False
        If True then drop `glon_peak` and `glat_peak` from table

    Returns
    -------
    bgps : pandas.DataFrame
        Output catalog
    """
    # Read in catalog
    bgps = _pd.read_csv(d.cat_dir + d.bgps_ext_filen.format('2.0.1', 'vel',
                        'csv'), na_values=['-1000.0'])
    # Column names for each molecule
    cols = ['vlsr_{0}_f'.format(s) for s in ('hco','nnh','cs','nh3')]
    for col in cols:
        bgps[col] = 0
    # Factor binary flag into four distinct 1/0 flags
    mask = _np.logical_not(bgps.vlsr_f.isin([-1,0]))
    bgps['vlsr_hco_f'][mask] = \
        bgps['vlsr_f'][mask].apply(lambda x: 1 if 1 in bin_factor(x) else 0)
    bgps['vlsr_nnh_f'][mask] = \
        bgps['vlsr_f'][mask].apply(lambda x: 1 if 2 in bin_factor(x) else 0)
    bgps['vlsr_cs_f'][mask] = \
        bgps['vlsr_f'][mask].apply(lambda x: 1 if 4 in bin_factor(x) else 0)
    bgps['vlsr_nh3_f'][mask] = \
        bgps['vlsr_f'][mask].apply(lambda x: 1 if 8 in bin_factor(x) else 0)
    # Collect velocities in one category
    bgps['all_vlsr'] = bgps['vlsr']
    grs_mask = (bgps.vlsr_f == 0) & (bgps.grs_vlsr_f == 1)
    bgps['all_vlsr'][grs_mask] = bgps['grs_vlsr'][grs_mask]
    # Don't clobber coordinates
    if no_clobber:
        bgps = bgps.drop(labels=['glon_peak', 'glat_peak'], axis=1)
    return bgps


def read_molcat():
    """
    Read BGPS HHT molecular line catalog. Citation: Shirley et al. (2013).

    Returns
    -------
    molcat : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    na_values = ['', '99.99', '-200.0', '-200.00', '-200.000', '999.9',
        '999.99', '999.999']
    molcat = _pd.read_csv(d.cat_dir + d.molcat_filen.format('main'),
        na_values=na_values)
    molcat = molcat.drop(labels=['mol_ir1', 'mol_ir2', 'mol_ir3', 'mol_ir4',
        'mol_c_l', 'mol_c_b'], axis=1)
    vlsr = _pd.read_csv(d.cat_dir + d.molcat_filen.format('vlsr'),
        na_values=na_values)
    vlsr = vlsr.drop(labels=['ref_f'], axis=1)
    gauss = _pd.read_csv(d.cat_dir + d.molcat_filen.format('gauss'),
        na_values=na_values)
    gauss = gauss.drop(labels=['hco_chisq', 'nnh_chisq'], axis=1)
    molcat = _pd.merge(molcat, vlsr, how='outer', on='cnum')
    molcat = _pd.merge(molcat, gauss, how='outer', on='cnum')
    molcat = pd_eq2gal(molcat, ['hht_ra','hht_dec'], ['hht_glon', 'hht_glat'])
    return molcat


def read_gbt_nh3(ret_idl=False):
    """
    Read BGPS GBT NH3 observation and temperature fit catalog. Citation:
    Dunham et al. (2011), Rosolowsky et al. (in prep.).

    Parameters
    ----------
    ret_idl : Boolean, default False
        Return list of record arrays from IDL save file

    Returns
    -------
    gbt_nh3 : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    idl_list : list
        IDL save file products
    """
    if type(ret_idl) != bool:
        raise TypeError('ret_idl must be Boolean')
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
    if ret_idl:
        return gbt_nh3, idl_list
    else:
        return gbt_nh3


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
    from besl.coord import sexstr2dec
    msx = _pd.read_csv(d.cat_dir + d.msx_filen)
    msx['ra'] = msx['ra_str'].apply(sexstr2dec, hd='h')
    msx['dec'] = msx['dec_str'].apply(sexstr2dec, hd='d')
    msx = pd_eq2gal(msx, ['ra', 'dec'], ['glon', 'glat'])
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


def read_ego(out_df_list=False):
    """
    Read EGO catalog. Citation: Cyganowski et al. (2008), Chen et al. (2013).

    Parameters
    ----------
    out_df_list : boolean, default False
        Condition to return list of individual DataFrames

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
    # cyganowski egos
    for i in range(1,6):
        df_list.append(_pd.read_csv(d.cat_dir + d.ego_filen.format(i,
        'cyganowski'), sep=';', skipinitialspace=True, skiprows=skip_list[i-1]))
    # chen egos
    chen_egos = _pd.read_csv(d.cat_dir + d.ego_filen.format(1, 'chen_coords'))
    chen_egos['ra'] = (chen_egos['ra_h'] + chen_egos['ra_m'] / 60. +
        chen_egos['ra_s'] / 60.**2) * 360. / 24.
    chen_egos['dec'] = chen_egos['dec_d'] + chen_egos['dec_m'] / 60. + \
        chen_egos['dec_s'] / 60.**2
    chen_egos = pd_eq2gal(chen_egos, labels=['ra', 'dec'], new_labels=['_Glon',
        '_Glat'])
    chen_egos = chen_egos[['EGO', '_Glon', '_Glat']]
    df_list.append(chen_egos)
    # bundle up out dataframe
    ego = _pd.concat(df_list, ignore_index=True)
    if out_df_list:
        return ego, df_list
    else:
        return ego


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
    mmb['dec'] = mmb.dec_d + _np.sign(mmb.dec_d) * \
                 (mmb.dec_m / 60. + mmb.dec_s / (60. * 60.))
    mmb = pd_eq2gal(mmb, ['ra', 'dec'], ['glon', 'glat'])
    return mmb


def read_pesta05(cat_type='csv'):
    """
    Read aggregate methanol maser survey catalog. Citation: Pestalozzi (2005).

    Parameters
    ----------
    cat_type : str, default 'csv'
        Type of catalog file to read, valid options include 'csv' and 'fit'. If
        'csv' read in as a `pandas.DataFrame`, else if 'fit' read in as a
        `astropy.io.fit.Fits` object.

    Returns
    -------
    pesta : pandas.DataFrame, astropy.io.fits.Fits
    """
    valid_cats = ['csv', 'fit']
    if cat_type not in ['csv', 'fit']:
        raise Exception('`cat_type` must be one of {0}'.format(valid_cats))
    path = d.cat_dir + d.pesta_metho_filen.format(cat_type)
    if cat_type == 'csv':
        return _pd.read_csv(path)
    elif cat_type == 'fit':
        return fits.open(path)
    else:
        raise Exception('Unexpected exception')


def read_gbt_h2o():
    """
    Read BGPS GBT H2O maser survey catalog. Citation: Svoboda et al. (2013).

    Returns
    -------
    gbt_h2o : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    tmb2jy = 0.463
    gbt_h2o = _pd.read_csv(d.cat_dir + d.gbt_h2o_filen,
        na_values=['-999.000000', '-9'], skiprows=21)
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
    from besl.coord import sexstr2dec
    arcetri_ces = _pd.read_csv(d.cat_dir + d.arcetri_ces_filen, sep=';',
        skipinitialspace=True, skiprows=1)
    arcetri_ces['ra'] = arcetri_ces['ra_str'].apply(sexstr2dec, hd='h')
    arcetri_ces['dec'] = arcetri_ces['dec_str'].apply(sexstr2dec, hd='d')
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
    arcetri_val['h2o_f'][arcetri_val.snr > 5] = 1
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


def read_cornish(exten='all'):
    """
    Read CORNISH survey catalog. Citation: Purcell et al. (2013).

    Parameters
    ----------
    exten : str, default 'all'
        CORNISH catalog subset. Valid types: all, hii, uchii.

    Returns
    -------
    corn : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    if exten not in ['all', 'uchii', 'hii']:
        raise ValueError
    corn = _pd.read_csv(d.cat_dir + d.cornish_filen.format(exten))
    corn = corn.rename(columns={'l_deg': 'glon', 'b_deg': 'glat', 'RA_deg':
           'ra', 'Dec_deg': 'dec'})
    return corn


def read_dpdf(v=2):
    """
    Read Distance Probability Distribution Functions. Citation:
    Ellsworth-Bowers et al. (2013).

    Parameters
    ----------
    v : number {1, 2, 21}, default 2
        Version of BGPS to use

    Returns
    -------
    dpdf : list
        List of `astropy.io.fits.FitsHDU`
    """
    versions = [1, 2, 21]
    if v not in versions:
        raise ValueError('v = {0}. Must be in {1}.'.format(v, versions))
    dpdf = fits.open(d.cat_dir + d.dpdf_filen.format(v))
    return dpdf


def read_emaf_dist(v=2):
    """
    Read EMAF calculated distances and parameters. Citation: Ellsworth-Bowers
    et al. (2013).

    Parameters
    ----------
    v : number {1, 2, 21}, default 2
        Version of BGPS to use

    Returns
    -------
    emaf : pandas.DataFrame
    """
    versions = [1, 2, 21]
    if v not in versions:
        raise ValueError('v = {0}. Must be in {1}.'.format(v, versions))
    if v == 1:
        emaf = _pd.read_csv(d.cat_dir + d.emaf_filen.format(v), sep=';',
            skipinitialspace=True, skiprows=33, na_values=['---'])
        return emaf
    if v == 2:
        emaf = _pd.read_csv(d.cat_dir + d.emaf_filen.format(v), sep=';',
            skipinitialspace=True, skiprows=31, na_values=['---'])
        return emaf
    if v == 21:
        emaf = _pd.read_csv(d.cat_dir + d.emaf_filen.format(v))
        return emaf


def read_oh94_dust(model_type='mrn', modeln=0):
    """
    Read OH dust opacity models. Citation: Ossenkopf & Henning (1994).

    Parameters
    ----------
    model_type : str, default 'mrn'
        OH model type, valid inputs: mrn, thick, thin
    modeln : number, default 0
        OH model density number, valid inputs: 0, 5, 6, 7, 8

    Returns
    -------
    oh94_dust : sp.interp1d
        Scipy interpolation function
    """
    if model_type not in ['mrn', 'thick', 'thin']:
        raise ValueError('Invalid model type.')
    if modeln not in [0, 5, 6, 7, 8]:
        raise ValueError('Invalid model density number.')
    x, y = _np.loadtxt(d.cat_dir + d.oh94.filen.format(model_type,
        modeln)).T
    oh94_dust = interp1d(x, y, kind='linear')
    return oh94_dust


### Region match
def select_bgps_field(lon, lat, coord_type='eq', bool_out=False):
    """
    Select BGPS image field name from sky position in galactic longitude and
    latitude.

    Parameters
    ----------
    glon : number
        Galactic longitude in decimal degrees
    glat : number
        Galactic latitude in decimal degrees
    coord_type : str, default 'eq'
        Choose 'gal' for Galactic or 'eq' for Equatorial
    bool_out : boolean
        Return True if in BGPS

    Returns
    -------
    field : str
        Image field name
    """
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
    if ((((lon > 359.4265) & (lon < 360)) | ((lon < 0.4515) & (lon > 0))) &
        ((lat < 0.629857) & (lat > -0.629857))):
        return 'l000'
    field = bgps_bounds[(bgps_bounds.glon_min < lon) &
                        (bgps_bounds.glon_max > lon) &
                        (bgps_bounds.glat_min < lat) &
                        (bgps_bounds.glat_max > lat)].field
    if bool_out:
        if len(field) == 0:
            return False
        elif len(field) > 0:
            return True
    if len(field) == 0:
        return _np.nan
    else:
        return field.values[0]

def num_in_bounds(cat, fields, cat_labels=['_Glon', '_Glat'],
    field_labels=['field', 'glon_min', 'glon_max', 'glat_min', 'glat_max'],
    sum_fields=False):
    """
    Calculate the number of sources in the overlap region of a catalog given
    field names and coordinates for the maximum and minimum in longitude and
    latitude of each field. Note that coordinate values should be in decimal
    degrees with longitude 0 - 360.

    Parameters
    ----------
    cat : pandas.DataFrame
        Catalog with targets to match
    fields : pandas.DataFrame
        Catalog with field names and edge coordinates
    cat_labels : list
        List of column labels for longitude and latitude coordinates in `cat`
    field_labels : list
        List of column labels for field-id, min longitude, max longitude, min
        latitude, and max latitude
    sum_fields : bool, default False
        Return the sum of all matches in the field

    Returns
    -------
    in_bounds : dict
        Dictionary with keys as field-id names in dataframe `fields` to the
        dataframe indices of overlapping sources in `cat`. If `sum_fields` is
        set to True, then a number is returned for the sum of all matched
        sources in the fields.
    """
    in_bounds = {}
    for ii in fields.index:
        field_id = fields.ix[ii, field_labels[0]]
        lon_min = fields.ix[ii, field_labels[1]]
        lon_max = fields.ix[ii, field_labels[2]]
        lat_min = fields.ix[ii, field_labels[3]]
        lat_max = fields.ix[ii, field_labels[4]]
        lon_width = ang_diff(lon_max, lon_min)
        lon_cen = _np.mod((lon_max + lon_min) / 2., 360)
        match = cat[(_np.abs(ang_diff(cat[cat_labels[0]].values, lon_cen)) <
                     lon_width / 2.) &
                    (cat[cat_labels[1]] > lat_min) &
                    (cat[cat_labels[1]] < lat_max)].index
        in_bounds[field_id] = match
    if sum_fields:
        return _np.sum([a.shape[0] for a in in_bounds.values()])
    else:
        return in_bounds

def attab2df(filen, **kwargs):
    """
    Read a CDS table and convert it into a `pandas.DataFrame`. Uses the
    `astropy.io.ascii` table reader and writers. Arguments are passed on to
    `astropy.io.ascii.read`.

    Parameters
    ----------
    filen : str
        Filename including extension.

    Returns
    -------
    df : pandas.DataFrame
    """
    tab = ascii.read(filen, **kwargs)
    with TemporaryFile() as tmp:
        ascii.write(tab, tmp, Writer=ascii.Tab)
        tmp.seek(0)  # rewind to file head
        df = _pd.read_csv(tmp, sep='\t', na_values=['--'])
    return df


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
    coord_type : str ('eq', 'gal'), default 'eq'
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
    if bgps_ver not in [1, 2]:
        raise ValueError('bgps_ver = {}. Must be 1 or 2.'.format(v))
    from besl.coord import nearest_match_coords
    search_rad = _np.sqrt(2) * pix_size / 2.
    index_list = []
    # read in catalog
    bgps = read_bgps(v=bgps_ver)
    field = bgps[bgps.cnum == cnum]['field'].values[0]
    fcnum = bgps[bgps.cnum == cnum]['field_cnum'].values[0]
    # read in label mask
    bgps_label = fits.open(d.bgps_dir +
        'v2.0_ds2_{}_13pca_labelmask.fits'.format(field))
    bgps_wcs = wcs.WCS(bgps_label[0].header)
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


### DS9 regions
def create_point_region(lon, lat, text=[], out_filen='ds9', marker='circle',
        coord_type='fk5', color='green', offset=5):
    """
    Create a DS9 region file from a list of longitude and latitude coordinates.

    Parameters
    ----------
    lon : array-like
        Longitude in decimal degrees
    lat : array-like
        Latitude in decimal degrees
    text : array-like
        List of text labels to place at offset position
    out_filen : str, default 'ds9.reg'
        Filename of output DS9 region file.
    marker : str, default 'circle'
        Region marker type. Supported: {circle, box, diamond, cross, x, arrow,
        boxcircle}. If 'none' no markers are plotted.
    coord_type : str, default 'fk5'
        Celestial coordinate type. Supported: {image, linear, fk4, fk5,
        galactic, ecliptic, icrs, physical, amplifier, detector}.
    color : str, default 'green'
        Region color. Supported: {white, black, red, green, blue, cyan, magenta,
        yellow}.
    offset : number, default 5 arcsec
        Offset for text from the point center

    Returns
    -------
    out_file : regions file
    """
    if marker not in ['circle', 'box', 'diamond', 'cross', 'x', 'arrow',
        'boxcircle', 'none']:
        raise ValueError('Invalid marker string.')
    if coord_type not in ['image', 'linear', 'fk4', 'fk5', 'galactic',
        'ecliptic', 'icrs', 'physical', 'amplifier', 'detector']:
        raise ValueError('Invalid coordinate type.')
    if color not in ['white', 'black', 'red', 'green', 'blue', 'cyan',
        'magenta', 'yellow']:
        raise ValueError('Invalid color.')
    if len(text) != len(lon):
        text = [''] * len(lon)
    out_file = open(out_filen + '.reg', 'w')
    # header file
    out_file.write('global color={0} font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0\n'.format(color))
    # coordinate system
    out_file.write('{0}\n'.format(coord_type))
    for i in xrange(lon.shape[0]):
        out_file.write(
        '{marker} point {lon} {lat} # text={lcb}{txt}{rcb}\n'.format(
            marker=marker,
            lon=lon[i],
            lat=lat[i],
            txt=text[i],
            lcb='{',
            rcb='}'))
    # TODO add text only option implementation
    out_file.close()
    print '-- DS9 region file written to {}.reg'.format(out_filen)
    return

