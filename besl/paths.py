"""
============
System Paths
============

Module for relative and absolute paths for files and scripts.

"""

import os as _os


### Directories, paths, and environment variables
class Dirs(object):
    """
    Object to hold directories for interactive editing of paths.
    """
    def __init__(self):
        self.root_dir = '/mnt/eld_data/'
        self.cat_dir = self.root_dir + 'Catalogs/'
        self.collected = '_collected/'
        self.bgps_dir = self.root_dir + 'BGPS/Images/v2.0.0/'
        self.working_dir = _os.getcwd() + '/'
        self.out_dir = self.working_dir + 'matched_cats/'
        self.bgps_filen = 'bgps/bgps_v{}.csv'
        self.bgps_ext_filen = 'bgps/bgps_v{}_{}.{}'
        self.bgps_bounds_filen = 'bgps/bgps_v2.0_bounds.csv'
        self.molcat_filen = 'bgps/bgps_molcat_{}.csv'
        self.wise_filen = 'wise/wise_0-90.csv'
        self.msx_filen = 'red_msx/rms_msx_urquhart.csv'
        self.robit_filen = 'red_robitaille/red_robitaille.csv'
        self.ego_filen = 'ego/ego_{}_{}.csv'
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
        self.dpdf_filen = 'emaf/BGPS_V{}_dpdf_table.fits'
        self.emaf_filen = 'emaf/emaf_dist_V{}.csv'
        self.gbt_nh3_1_filen = 'bgps/nh3_11B48_fit_objects.sav'
        self.gbt_nh3_2_filen = 'bgps/nh3_bgps_fit_objects.sav'
        self.gbt_nh3_3_filen = 'bgps/nh3_rms_fit_objects.sav'
        self.oh94_filen = 'oh94_dust/{}{}.asc'
        self.cornish_filen = 'cornish/cornish_{}.csv'
        self.pesta_metho_filen = 'pestalozzi05/pestalozzi05.{}'


all_paths = Dirs()


