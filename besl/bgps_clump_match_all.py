"""
===============
Clump Match All
===============

Merge catalogs based on clump label masks

"""

import os
import numpy as _np
import pandas as _pd
import catalog
from .catalog import (read_bgps, read_cat)
from .image import sample_bgps_img


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
        if verbose:
            print '-- clump {:>4d}'.format(cnum)
        # match egos
        if (glat < 1.05) & (glat > -1.05) & (glon < 65):
            ego_match_list = catalog.clump_match(ego_hs, cnum,
                coord_type='gal')
            bgps['ego_n'][cnum_select] = len(ego_match_list)
        # match robit
        if (glat < 65):
            robit_agb_match_list = catalog.clump_match(robit_hs, cnum,
                coord_type='gal')
            robit_yso_match_list = catalog.clump_match(robit_hs, cnum,
                                                   coord_type='gal')
            bgps['robit_agb_n'][cnum_select] = len(robit_agb_match_list)
            bgps['robit_yso_n'][cnum_select] = len(robit_yso_match_list)
        # match rms msx
        if (glat < 5) & (glat > -5) & (glon > 10) & (glat < 220):
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


###############################################################################
#                            Class Based Approach
###############################################################################


class Matcher(object):
    """
    Matching class to match sources sources to the BGPS based on the label
    masks. A data object is given and then the catalog is processed.

    Parameters
    ----------
    The `data` object is required to have the following attributes:
    name : str
        Name used for a unique identifier for output
    cat : pd.DataFrame
        The catalog to match
    lon_col : str
        Galactic longitude column name in decimal degrees
    lat_col : str
        Galactic latitude column name in decimal degrees

    Attributes
    ----------
    """
    v = 210
    cnum_col = 'v210cnum'

    def __init__(self, data):
        # Parameters from data object
        self.name = data.name
        self.cat = data.cat
        self.lon_col = data.lon_col
        self.lat_col = data.lat_col
        self.det_col = data.det_col
        self.det_flags = data.det_flags
        self.choose_col = data.choose_col
        self.noise_col = data.noise_col
        self.data = data
        # BGPS data
        self.bgps = read_bgps(v=self.v).set_index(self.cnum_col)
        # Process and match
        self._add_new_cols()
        self._make_haystack()
        self._match()

    def _add_new_cols(self):
        """
        Add new columns to the BGPS catalog depending on the available
        columns to match to in the child catalog.
        """
        # Don't clobber original columns in BGPS
        #  Do first so as to not rename following named flags
        for col in self.cat.columns:
            if col in self.bgps.columns:
                self.bgps['_' + col] = _np.nan
            else:
                self.bgps[col] = _np.nan
        # New column for number of matched sources
        self.bgps_count_col = self.name + '_n'
        self.bgps[self.bgps_count_col] = _np.nan
        # For number of detections
        if self.det_col is not None:
            self.bgps_det_col = self.name + '_f'
            self.bgps[self.bgps_det_col] = _np.nan

    def _make_haystack(self):
        """
        Make 2xN array of coordinate positions for each source in the child
        catalog.
        """
        self.haystack = self.cat[[self.lon_col,
                                  self.lat_col]].values

    def _match(self):
        """
        Match each source in the child catalog to the BGPS.
        """
        self.matched_ix = {}
        ids = self.cat.index
        haystack = self.haystack
        for cat_ix, coord in zip(ids, haystack):
            cnum = sample_bgps_img(coord[0], coord[1], v=self.v)
            self._enter_matched(cat_ix, cnum)

    def _enter_matched(self, ix, cnum):
        """
        Insert the indices of the matched sources into the matched source
        dictionary with the BGPS cnum as the key.
        """
        if (not _np.isnan(cnum)) & (cnum != 0):
            self.matched_ix.setdefault(cnum, [])
            self.matched_ix[cnum].append(ix)

    def process(self):
        """
        Simple processing to add:
            number of matches
            number of detections
                If `det_col` is specified.
            all columns
                If `choose_col` and `noise_col` are specified then the
                maximum in `choose_col` will be used, or if all are null
                values, then source with minimum noise will be chosen.
                If `noise_col` is `None` then the first source in the sub-index
                will be used.
        columns to BGPS catalog.
        """
        # New column for number of matched sources
        for cnum, cat_indices in self.matched_ix.iteritems():
            self.bgps.ix[cnum, self.bgps_count_col] = len(cat_indices)
            if self.det_col is not None:
                matches = self.cat.ix[cat_indices, self.det_col]
                num_dets = matches[matches.isin(self.det_flags)].shape[0]
                self.bgps.ix[cnum, self.bgps_det_col] = num_dets
            if self.choose_col is not None:
                choose_ix = self.cat.ix[cat_indices, self.choose_col].idxmax()
                if _np.isnan(choose_ix) & (self.noise_col is not None):
                    choose_ix = self.cat.ix[cat_indices,
                                            self.noise_col].idxmin()
                else:
                    choose_ix = cat_indices[0]
                self.bgps.ix[cnum, self.cat.columns] = self.cat.ix[choose_ix]

    def to_csv(self):
        """
        Write BGPS catalog to `.csv` file.
        """
        self.bgps = self.bgps.rename(columns={col: self.name + '_' + col for
                                              col in self.bgps.columns[20:]})
        self.bgps.to_csv('bgps_' + self.name + '.csv', index=False)


class DataSet(object):
    def __init__(self):
        self.all_data = []
        all_objs = [WaterGbt,
                    WaterArcetri,
                    WaterHops,
                    WaterRms,
                    Cornish,
                    Egos,
                    AmmoniaGbt,
                    MethoPandian,
                    MethoPestalozzi,
                    MethoMmb,
                    Higal70,
                    RedSpitzer,
                    RedMsx,
                    Molcat]
        for obj in all_objs:
            data = obj()
            data.match()
            data.write()
            self.all_data.append(data)
        self.process()

    def _merge(self):
        print '-- Merging data'
        merged_data = read_bgps(v=210).set_index('v210cnum')
        for data in self.all_data:
            merged_data = merged_data.merge(data, how='left',
                                            left_index=True,
                                            right_index=True)
        self.merged_data = merged_data

    def _append_evo_flags(self):
        print '-- Adding evolutionary flags'
        self.merged_data = append_evo_flags(bgps=self.merged_data)

    def _write(self):
        print '-- Writing all merged data'
        self.merged_data.to_csv('bgps_v210_all.csv', index=False)

    def process(self):
        self._merge()
        self._append_evo_flags()
        self._write()


class Data(object):
    """
    Parent class for object-catalogs to be matched with `Matcher`
    """
    def match(self):
        print '-- Matching {0}'.format(self.name)
        self.matcher = Matcher(self)
        self.matcher.process()

    def write(self):
        print '-- Writing {0}'.format(self.name)
        self.matcher.to_csv()


def append_evo_flags(bgps):
    """
    Calculate and append evolutionary flags to BGPS catalog

    Parameters
    ----------

    Returns
    -------
    """
    evo_flags = ['h2o_f', 'ch3oh_f', 'ego_f', 'ir_f', 'uchii_f', 'sf_f']
    for col in evo_flags:
        bgps[col] = _np.nan
    # H2O flags
    bgps[((bgps['h2o_gbt_n'] > 0) & (bgps['h2o_gbt_f'] == 0)) &
         _np.logical_not(bgps['h2o_arc_f'] > 0) &
         _np.logical_not(bgps['h2o_rms_n'] > 0) &
         _np.logical_not(bgps['h2o_hops_n'] > 0)]['h2o_f'] = 0
    bgps[(bgps['h2o_gbt_f'] > 0) |
         (bgps['h2o_rms_n'] > 0) |
         (bgps['h2o_arc_n'] > 0) |
         (bgps['h2o_hops_n'] > 0)]['h2o_f'] = 1
    # CH3OH flags
    bgps[(bgps['ch3oh_pesta'] > 0) |
         (bgps['ch3oh_pandi'] > 0) |
         (bgps['ch3oh_mmb'] > 0)]['ch3oh_f'] = 1
    # EGO flags
    bgps[(bgps['ego_n'] > 0)]['ego_f'] = 1
    # IR flags
    bgps[(bgps['robit_f'] > 0) |
         (bgps['red_msx_f'] > 0) |
         (bgps['ego_n'] > 0)]['ir_f'] = 1
    bgps[(bgps['ir_f'] != 1) &
         (((bgps['robit_n'] > 0) & (bgps['robit_f'] == 0)) |
          ((bgps['red_msx_n'] > 0) & (bgps['red_msx_f'] == 0)))]['ir_f'] = 2
    # UCHII flags
    bgps[(bgps['corn_n'] > 0)]['uchii_f'] = 1
    # Starless
    bgps[(bgps['h2o_f'] == 0) &
         (bgps['ch3oh_f'] != 1) &
         _np.logical_not(bgps['ir_f'].isin([1, 2])) &
         (bgps['uchii_f'] != 1)]['sf_f'] = 0
    bgps[(bgps['h2o_f'] == 1) |
         (bgps['ch3oh_f'] == 1) |
         (bgps['ir_f'] == 1) |
         (bgps['uchii_f'] == 1)]['sf_f'] = 1
    return bgps


###############################################################################
#                          Catalog Data Objects
###############################################################################


class WaterGbt(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'h2o_gbt'
        self.cat = read_cat('gbt_h2o')
        self.lon_col = 'h2o_glon'
        self.lat_col = 'h2o_glat'
        self.det_col = 'h2o_f'
        self.det_flags = [1]
        self.choose_col = 'h2o_tpk'
        self.noise_col = 'h2o_tpk_err'


class WaterArcetri(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'h2o_arc'
        self.cat = read_cat('valdettaro01_arcetri')
        self.lon_col = '_Glon'
        self.lat_col = '_Glat'
        self.det_col = 'h2o_f'
        self.det_flags = [1]
        self.choose_col = 'Stot'
        self.noise_col = 'Sig'


class WaterHops(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'h2o_hops'
        self.cat = read_cat('walsh11_hops_h2o')
        self.lon_col = 'lPeak_deg'
        self.lat_col = 'bPeak_deg'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 'T_peak_K'
        self.noise_col = 'RMS_K'


class WaterRms(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'h2o_rms'
        self.cat = read_cat('urquhart11_red_msx_h2o')
        self.lon_col = '_Glon_1_'
        self.lat_col = '_Glat_1_'
        self.det_col = 'H2O_1_'
        self.det_flags = ['y']
        self.choose_col = 'log_SdV__2_'
        self.noise_col = 'rms_2_'


class Cornish(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'corn'
        self.cat = read_cat('cornish_uchii')
        self.lon_col = 'l_deg'
        self.lat_col = 'b_deg'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 'Flux_mJy'
        self.noise_col = 'dFlux_mJy'


class Egos(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'ego'
        self.cat = read_cat('ego_all')
        self.lon_col = '_Glon'
        self.lat_col = '_Glat'
        self.det_col = None
        self.det_flags = None
        self.choose_col = '[4.5]'
        self.noise_col = None


class AmmoniaGbt(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'nh3_gbt'
        self.cat = read_cat('gbt_nh3')
        self.lon_col = 'glon'
        self.lat_col = 'glat'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 'pk11'
        self.noise_col = 'noise11'


class MethoPandian(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'ch3oh_pandi'
        self.cat = read_cat('pandian11')
        self.lon_col = 'glon'
        self.lat_col = 'glat'
        self.det_col = None
        self.det_flags = None
        self.choose_col = None
        self.noise_col = None


class MethoPestalozzi(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'ch3oh_pesta'
        self.cat = read_cat('pestalozzi05')
        self.lon_col = '_Glon'
        self.lat_col = '_Glat'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 'PFlux'
        self.noise_col = None


class MethoMmb(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'ch3oh_mmb'
        self.cat = read_cat('mmb_all')
        self.lon_col = 'glon'
        self.lat_col = 'glat'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 'spk_mx'
        self.noise_col = None


class Higal70(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'higal70'
        self.cat = read_cat('higal_70')
        self.lon_col = 'glon_70'
        self.lat_col = 'glat_70'
        self.det_col = None
        self.det_flags = None
        self.choose_col = 's70cx'
        self.noise_col = 'err_s70cx'


class RedSpitzer(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'robit'
        self.cat = read_cat('robitaille08_red_spitzer')
        self.lon_col = '_Glon'
        self.lat_col = '_Glat'
        self.det_col = 'Class'
        self.det_flags = ['cYSO']
        self.choose_col = '[8.0]'
        self.noise_col = None


class RedMsx(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'red_msx'
        self.cat = read_cat('lumsden13_red_msx')
        self.lon_col = 'glon'
        self.lat_col = 'glat'
        self.det_col = 'Type'
        self.det_flags = ['YSO', 'HII/YSO', 'Young/old star']
        self.choose_col = None
        self.noise_col = None


class Molcat(Data):
    def __init__(self):
        # Catalog parameters
        self.name = 'mol'
        self.cat = read_cat('shirley13_molcat')
        self.lon_col = 'hht_glon'
        self.lat_col = 'hht_glat'
        self.det_col = 'mol_vlsr_f'
        self.det_flags = [1, 2, 3]
        self.choose_col = 'hco_tpk'
        self.noise_col = 'hco_tpk_err'


