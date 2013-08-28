"""
===============
Clump Match All
===============

Merge catalogs based on clump label masks

"""

import os
import catalog
import numpy as _np
import pandas as _pd
import ipdb as pdb

def clump_match_water(bgps=[], out_filen='bgps_maser', verbose=False):
    """
    Match maser catalog observations to the BGPS. Includes BGPS GBT, Red MSX,
    Arcetri, MMB, and HOPS.

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_maser'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    gbt_h2o = catalog.read_gbt_h2o()
    rms_h2o = catalog.read_rms_h2o()
    arc_val = catalog.read_arcetri_valdettaro()
    hops = catalog.read_hops()
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # add new columns
    new_cols = ['h2o_gbt_f', 'h2o_gbt_n', 'h2o_arc_f', 'h2o_arc_n',
            'h2o_hops_f', 'h2o_rms_f', 'h2o_rms_n']
    gbt_cols = gbt_h2o.columns.drop(labels=['h2o_glon', 'h2o_glat', 'h2o_f'])
    for col in new_cols:
        bgps[col] = _np.nan
    for col in gbt_cols:
        bgps[col] = _np.nan
    # make haystacks
    gbt_h2o_hs = gbt_h2o[['h2o_glon', 'h2o_glat']].values # galactic
    rms_h2o_hs = rms_h2o[['_Glon_y', '_Glat_y']].values # galactic
    arc_val_hs = arc_val[['_Glon', '_Glat']].values # galactic
    hops_hs = hops[['lWeight_deg', 'bWeight_deg']].values # galactic
    # loop through clumps
    for cnum in bgps['cnum']:
        cnum_select = bgps.cnum == cnum
        c_index = _np.argwhere(cnum_select)[0][0]
        glat = bgps[cnum_select].glat_cen.values[0]
        glon = bgps[cnum_select].glon_cen.values[0]
        c_ra = bgps[cnum_select].ra.values[0]
        c_dec = bgps[cnum_select].dec.values[0]
        # match hops
        if ((glat < 0.5) & (glat > -0.5) & ((glon > 290) | (glon < 30))):
            hop_match_list = catalog.clump_match(hops_hs, cnum,
                coord_type='gal')
            bgps['h2o_hops_f'][cnum_select] = len(hop_match_list)
        # match bgps gbt
        gbt_match_list = catalog.clump_match(gbt_h2o_hs, cnum,
            coord_type='gal')
        h2o_gbt_num_detects = _np.sum(gbt_h2o.h2o_f.ix[gbt_match_list])
        bgps['h2o_gbt_n'][cnum_select] = len(gbt_match_list)
        bgps['h2o_gbt_f'][cnum_select] = h2o_gbt_num_detects
        if h2o_gbt_num_detects > 0:
            max_index = gbt_h2o['h2o_tpk'].ix[gbt_match_list].argmax()
            bgps.ix[c_index, gbt_cols] = \
                gbt_h2o.ix[gbt_match_list[max_index]]
        # match rms h2o
        rms_match_list = catalog.clump_match(rms_h2o_hs, cnum,
            coord_type='gal')
        bgps['h2o_rms_n'][cnum_select] = len(rms_match_list)
        bgps['h2o_rms_f'][cnum_select] = \
            _np.sum(rms_h2o.h2o_f.ix[rms_match_list])
        # match arcetri
        arc_match_list = catalog.clump_match(arc_val_hs, cnum,
            coord_type='gal')
        bgps['h2o_arc_n'][cnum_select] = len(arc_match_list)
        bgps['h2o_arc_f'][cnum_select] = \
            _np.sum(arc_val.h2o_f.ix[arc_match_list])
        if verbose:
            print '-- clump {:>4d}'.format(cnum)
    bgps['h2o_f'] = _np.nan
    bgps['h2o_f'][(bgps.h2o_gbt_f > 0) | (bgps.h2o_arc_f > 0) |
        (bgps.h2o_rms_f > 0) | (bgps.h2o_hops_f > 0)] = 1
    bgps['h2o_f'][(bgps.h2o_f != 1) & ((bgps.h2o_gbt_f == 0) & (bgps.h2o_gbt_n >
        0))] = 0
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- Maser catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_hii(bgps=[], out_filen='bgps_hii', verbose=False):
    """
    Match HII and UCHII catalog observations to the BGPS. Include CORNISH and
    HRDS.

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_hii'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    corn = catalog.read_cornish(exten='hii')
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # add new columns
    new_cols = ['corn_n']
    for col in new_cols:
        bgps[col] = _np.nan
    # make haystacks
    corn_hs = corn[['glon', 'glat']].values # galactic
    # loop through clumps
    for cnum in bgps['cnum']:
        cnum_select = bgps.cnum == cnum
        glat = bgps[cnum_select].glat_cen.values[0]
        glon = bgps[cnum_select].glon_cen.values[0]
        c_ra = bgps[cnum_select].ra.values[0]
        c_dec = bgps[cnum_select].dec.values[0]
        if verbose:
            print '-- clump {:>4d}'.format(cnum)
        # match cornish
        if (glat < 1.0) & (glat > -1.0) & (glon > 9.95) & (glon < 65.55):
            corn_match_list = catalog.clump_match(corn_hs, cnum,
                coord_type='gal')
            bgps['corn_n'][cnum_select] = len(corn_match_list)
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- Hii catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_ir(bgps=[], out_filen='bgps_ir', verbose=False):
    """
    Match IR point source catalog observations to the BGPS. Includes EGO, RMS,
    and Robitaille.

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_ir'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    ego = catalog.read_ego()
    robit = catalog.read_robitaille()
    msx = catalog.read_msx()
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # add new columns
    new_cols = ['ego_n', 'msx_n', 'robit_n']
    for col in new_cols:
        bgps[col] = _np.nan
    # make haystacks
    ego_hs = ego[['_Glon', '_Glat']].values # galactic
    robit_hs = robit[['_Glon', '_Glat']].values # galactic
    msx_hs = msx[['ra', 'dec']].values # equatorial
    # loop through clumps
    for cnum in bgps['cnum']:
        cnum_select = bgps.cnum == cnum
        glat = bgps[cnum_select].glat_cen.values[0]
        glon = bgps[cnum_select].glon_cen.values[0]
        c_ra = bgps[cnum_select].ra.values[0]
        c_dec = bgps[cnum_select].dec.values[0]
        if verbose:
            print '-- clump {:>4d}'.format(cnum)
        # match egos
        if ((glat < 1.05) & (glat > -1.05) & (glon > 10) & (glat < 65)):
            ego_match_list = catalog.clump_match(ego_hs, cnum,
                coord_type='gal')
            bgps['ego_n'][cnum_select] = len(ego_match_list)
        # match robit
        if ((glat < 1.05) & (glat > -1.05) & (glon > 10) & (glat < 65)):
            robit_match_list = catalog.clump_match(robit_hs, cnum,
                coord_type='gal')
            bgps['robit_n'][cnum_select] = len(robit_match_list)
        # match rms msx
        if ((glat < 5) & (glat > -5) & (glon > 10) & (glat < 220)):
            msx_match_list = catalog.clump_match(msx_hs, cnum,
                coord_type='eq')
            bgps['msx_n'][cnum_select] = len(msx_match_list)
        # TODO add red wise
    bgps['ir_f'] = -9
    bgps['ir_f'][(bgps.ego_n > 0) | (bgps.msx_n > 0) |
        (bgps.robit_n > 0)] = 1
    bgps['ir_f'][(bgps.ego_n == 0) & (bgps.msx_n == 0) &
        (bgps.robit_n == 0)] = 0
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- IR catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_molcat(bgps=[], out_filen='bgps_molcat', verbose=False):
    """
    Match the BGPS HCO+/N2H+ molecular line survey observations to the BGPS.
    Citation: Shirley et al. (2013).

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_molcat.csv'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    pdb.set_trace()
    molcat = catalog.read_molcat()
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # add new columns, molcat cnum clobbers bgps cnum
    molcat = molcat.rename(columns={'cnum': 'v1cnum'})
    mol_type = {'HCOP': 'hco_tpk', 'N2HP': 'nnh_tpk'}
    # add column for intensity of dominant molecule
    molcat['mol_int'] = _np.nan
    for i in molcat['mol_vlsr_src'].index:
        mol_select = molcat.ix[i, 'mol_vlsr_src']
        if mol_select in mol_type.keys():
            mol_int = molcat.ix[i, mol_type[mol_select]]
            molcat.ix[i, 'mol_int'] = mol_int
    # columns
    molcat_cols = molcat.columns
    for col in molcat_cols:
        bgps[col] = _np.nan
    bgps['mol_mult_n'] = _np.nan
    bgps['mol_mult_f'] = _np.nan
    bgps['mol_mult_vsep'] = _np.nan
    # make haystacks
    molcat_hs = molcat[['hht_ra', 'hht_dec']].values # galactic
    # loop through clumps
    for cnum in bgps['cnum']:
        cnum_select = bgps.cnum == cnum
        c_index = _np.argwhere(cnum_select)[0][0]
        glat = bgps.ix[c_index, 'glat_cen']
        glon = bgps.ix[c_index, 'glon_cen']
        # match molcat
        # TODO mark -9 in outer regions which not observed in v1
        if (glon > 7.5) & (glon < 195):
            molcat_match_list = catalog.clump_match(molcat_hs, cnum,
                coord_type='eq')
            bgps['mol_mult_n'][cnum_select] = len(molcat_match_list)
            if verbose:
                print '-- {} matched to {} pointings'.format(cnum,
                    len(molcat_match_list))
            if len(molcat_match_list) == 1:
                bgps.ix[c_index, molcat_cols] = molcat.ix[molcat_match_list[0]]
                bgps['mol_mult_f'][cnum_select] = 0
                if molcat.ix[molcat_match_list]['mol_vlsr_f'] == 2:
                    bgps['mol_mult_f'][cnum_select] = 1
            elif len(molcat_match_list) > 1:
                flags = molcat.ix[molcat_match_list]['mol_vlsr_f'].values
                # if multiple component in a single spectra then confused
                if 2 in flags:
                    bgps['mol_mult_f'][cnum_select] = 1
                # if only single detection then not confused
                elif flags[(flags == 1) | (flags == 3)].shape[0] <= 1:
                    bgps['mol_mult_f'][cnum_select] = 0
                # if only non-detections
                elif flags[flags == 0].shape[0] == len(molcat_match_list):
                    bgps['mol_mult_f'][cnum_select] = 0
                    print molcat[molcat.mol_vlsr_f ==
                        0].ix[molcat_match_list]['mol_vlsr_f']
                # if more than one detection or self-absorbed and the
                # the velocity seperation between two components is more than
                # 5 km/s mark as confused, otherwise, not confused
                elif flags[(flags == 1) | (flags == 3)].shape[0] > 1:
                    vmin = molcat[(molcat.mol_vlsr_f == 1) |
                        (molcat.mol_vlsr_f ==
                            3)].ix[molcat_match_list]['mol_vlsr'].min()
                    vmax = molcat[(molcat.mol_vlsr_f == 1) |
                        (molcat.mol_vlsr_f == 3)].ix[molcat_match_list]['mol_vlsr'].max()
                    vsep = _np.abs(vmin - vmax)
                    bgps['mol_mult_vsep'][cnum_select] = vsep
                    if vsep < 5:
                        bgps['mol_mult_f'][cnum_select] = 0
                    else:
                        bgps['mol_mult_f'][cnum_select] = 1
                else:
                    raise Exception("Unexpected number of flags")
                # match values for component with peak intensity
                max_index = molcat['mol_int'].ix[molcat_match_list].argmax()
                bgps.ix[c_index, molcat_cols] = molcat.ix[molcat_match_list[max_index]]
                if verbose:
                    print bgps.ix[c_index, molcat_cols]
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- Molcat catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_gbt_nh3(bgps=[], out_filen='bgps_nh3', verbose=False):
    """
    Match the BGPS GBT NH3 survey observations to the BGPS.  Citation:
    Dunham et al. (2011), Rosolowsky et al. (in prep.).

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_nh3.csv'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    nh3 = catalog.read_gbt_nh3()
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # add new columns, molcat cnum clobbers bgps cnum
    nh3_cols = ['GLON', 'GLAT', 'TKIN', 'TKIN_ERR', 'VLSR', 'PK11', 'NOISE11',
                'PK22', 'NOISE22', 'PK33', 'NOISE33']
    nh3 = nh3[nh3_cols]
    nh3 = nh3.rename(columns={s: 'nh3_' + s.lower() for s in nh3_cols})
    nh3_cols = ['nh3_' + s.lower() for s in nh3_cols]
    for col in nh3_cols:
        bgps[col] = _np.nan
    bgps['nh3_mult_n'] = _np.nan
    # make haystacks
    nh3_hs = nh3[['nh3_glon', 'nh3_glat']].values # galactic
    # loop through clumps
    for cnum in bgps['cnum']:
        cnum_select = bgps.cnum == cnum
        c_index = _np.argwhere(cnum_select)[0][0]
        glat = bgps[cnum_select].glat_cen.values[0]
        glon = bgps[cnum_select].glon_cen.values[0]
        c_ra = bgps[cnum_select].ra.values[0]
        c_dec = bgps[cnum_select].dec.values[0]
        # match gbt nh3
        if verbose:
            print '-- clump {:>4d}'.format(cnum)
        if (glon > 7.5) & (glon < 200):
            nh3_match_list = catalog.clump_match(nh3_hs, cnum,
                coord_type='gal')
            bgps['nh3_mult_n'][cnum_select] = len(nh3_match_list)
            if len(nh3_match_list) > 0:
                max_index = nh3['nh3_pk11'].ix[nh3_match_list].argmax()
                bgps.ix[c_index, nh3_cols] = \
                    nh3.ix[nh3_match_list[max_index]]
                # TODO check for multiple components
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- NH3 Catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_metho(bgps=[], out_filen='bgps_metho', verbose=False):
    """
    Match known CH3OH maser catalogs to the BGPS.  Citation: Pandian et al.
    (2007, 2011), Pestalozzi et al. (2005), Caswell et al. (2010), Green et al.
    (2010).

    Paramters
    ---------
    bgps : pandas.DataFrame, default []
        BGPS catalog to match to, defaults to read vanilla catalog
    out_filen : string, default 'bgps_nh3.csv'
        Name of output catalog, comma seperated
    verbose : boolean, default False
        Print clump and number of matches

    Returns
    -------
    bgps : pd.DataFrame
    """
    # read in catalogs
    pandi = catalog.read_cat('pandian11')
    pesta = catalog.read_cat('pestalozzi05')
    mmb = catalog.read_mmb()
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # use general match TODO put BGPS label masks in cache
    bgps = clump_match_gen(pandi, bgps=bgps, prefix='pandi',
            verbose=verbose)
    bgps = clump_match_gen(pesta, bgps=bgps, prefix='pesta',
            verbose=verbose)
    bgps = clump_match_gen(mmb, bgps=bgps, coord_labels=['glon', 'glat'],
            prefix='mmb', verbose=verbose)
    # mark master ch3oh flag
    bgps['ch3oh_f'] = _np.nan
    bgps['ch3oh_f'][(bgps.pandi_n > 0) | (bgps.pesta_n > 0) |
                    (bgps.mmb_n > 0)] = 1
    # print to file
    bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
    print '-- Catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_gen(cat, bgps=[], coord_labels=['_Glon', '_Glat'],
    prefix='', out_filen=None, coord_type='gal', verbose=False,
    det_col=None, best_col=None, add_all_cols=False):
    """
    Match the BGPS to a general catalog with coordinate columns _Glon and
    _Glat.

    Parameters
    ----------
    cat : pd.DataFrame
        Catalog to match
    bgps : pd.DataFrame, default []
        BGPS or matched catalog to merge into. If empty list is passed then the
        default BGPS v2 is read in.
    prefix : str, default ''
        Name to append to beginning of columns from cat
    out_filen : str, default None
        Name of output file, if left None, no output file is written
    coord_type : str, default 'gal'
        Coordinate type, either Galactic 'gal' or Equatorial 'eq'
    verbose : True
        Print progress per clump
    det_col : str, default None
        Column name for detection flags, valid for 0 or 1. If `None` then do
        nothing.
    best_col : str, default None
        Column to use max value from when discriminating multiple matches if
        `add_all_cols` is set to True.
    add_all_cols : bool, default False
        Join all columns from input catalog `cat` to the BGPS.

    Returns
    -------
    bgps : pd.DataFrame
        Merged catalog
    """
    if cat[coord_labels[0]].min() < 0:
        raise ValueError('Longitude must be from 0 to 360')
    if len(bgps) == 0:
        bgps = catalog.read_bgps()
    # rename columns
    for col in cat.columns:
        cat = cat.rename(columns={col: prefix + '_' + col})
    # make sure not to clobber column names
    if _np.any(_np.in1d(cat.columns, bgps.columns)):
        overlap_cols = cat.columns[_np.in1d(cat.columns, bgps.columns)]
        for col in overlap_cols:
            cat = cat.rename(columns={col: '_' + col})
    # assign new columns to empty values
    for col in cat.columns:
        bgps[col] = _np.nan
    bgps[prefix + '_n'] = _np.nan
    if det_col is not None:
        bgps[prefix + '_f'] = _np.nan
    if add_all_cols:
        for col in cat.columns:
            bgps[col] = _np.nan
    # make haystack
    coord_labels = [prefix + '_' + col for col in coord_labels]
    cat_hs = cat[coord_labels].values
    # loop through clumps
    for cnum in bgps['cnum']:
        # BGPS properties
        cnum_select = bgps.cnum == cnum
        c_index = _np.argwhere(cnum_select)[0][0]
        # match cat
        match_list = catalog.clump_match(cat_hs, cnum, coord_type=coord_type)
        bgps[prefix + '_n'][cnum_select] = len(match_list)
        if det_col is not None:
            bgps[prefix + '_f'][cnum_select] = \
                cat[det_col].ix[match_list].sum()
        if add_all_cols:
            max_index = cat[best_col].ix[match_list].argmax()
            bgps.ix[c_index, cat.columns] = cat.ix[match_list[max_index]]
        if verbose:
            print '-- clump {0:>4d} : {1:>4d}'.format(cnum, len(match_list))
    if out_filen is not None:
        bgps.to_csv(os.getcwd() + '/' + out_filen + '.csv', index=False)
        print '-- Catalog file written to {}.csv'.format(out_filen)
    return bgps

def clump_match_all():
    """
    Match BGPS to all evolutionary indicators. Matches to NH3, HCO+/N2H+,
    H2O/CH3OH, IR, and HII.

    Returns
    -------
    bgps : pd.DataFrame
    """
    bgps = catalog.read_bgps()
    bgps_all = bgps.copy()
    df_list = []
    fn_list = [clump_match_molcat,
               clump_match_gbt_nh3,
               clump_match_water,
               clump_match_metho,
               clump_match_ir,
               clump_match_hii]
    for fn in fn_list:
        df_list.append(fn())
    for df in df_list:
        bgps_all = _pd.merge(bgps_all, df, how='outer')
    bgps_all.to_csv('bgps_all.csv', index=False)
    bgps_all.save('bgps_all.pickle')
    return bgps_all


