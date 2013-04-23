"""
Extract peak values from BGPS maps
"""

import numpy as _np
import pyfits as _pyfits
import pywcs as _pywcs

class Dirs(object):
    """
    Object to hold directories for interactive path editing
    """
    def __init__(self):
        self.root_dir = '/mnt/eld_data/'
        self.bgps1_img_dir = self.root_dir + 'BGPS/Images/v1.0.2/all/'
        self.bgps2_img_dir = self.root_dir + 'BGPS/Images/v2.0.0/'
        self.bgps1_img_filen = 'v1.0.2_{}_13pca_map50.fits'
        self.bgps1_med_filen = 'v1.0.2_{}_13pca_medmap50.fits'
        self.bgps1_rms_filen = 'v1.0.2_{}_13pca_noisemap50.fits'
        self.bgps2_img_filen = 'v2.0_ds2_{}_13pca_map20.fits'
        self.bgps2_med_filen = 'v2.0_ds2_{}_13pca_medmap20.fits'
        self.bgps2_rms_filen = 'v2.0_ds2_{}_13pca_noisemap20.fits'
d = Dirs()

def extract_peak_bgps_props(out_filen='bgps_pk_extract'):
    """
    Extract peak flux and noise values from the BGPS maps in Jy/beam. Citation:
    Ginsburg et al. (2013).

    Parameters
    ----------
    out_filen : string
        output catalog file in CSV format

    Returns
    -------
    molcat_pk : pandas.DataFrame
        Output catalog in a pandas DataFrame object
    """
    # TODO add extraction from v1 maps
    from besl.coord import eq2gal
    from besl.catalog import read_bgps, read_bgps_bounds, read_molcat, \
                             select_bgps_field
    # read in catalogs
    bgps_bounds = read_bgps_bounds()
    molcat = read_molcat()
    # add new columns
    new_cols = ['map_pk', 'medmap_pk', 'map_rms']
    for col in new_cols:
        molcat[col] = _np.nan
    # match clumps
    for clump in molcat.cnum.values:
        # index parameters
        cindex = _np.argwhere(molcat.cnum == clump)[0][0]
        ra, dec = molcat.ix[cindex]['hht_ra'], molcat.ix[cindex]['hht_dec']
        glon, glat = eq2gal(ra, dec)
        # bgps v2
        field = select_bgps_field(ra, dec, coord_type='eq')
        bgps2_img = _pyfits.open(d.bgps2_img_dir + d.bgps2_img_filen.format(field))
        bgps2_med = _pyfits.open(d.bgps2_img_dir + d.bgps2_med_filen.format(field))
        bgps2_rms = _pyfits.open(d.bgps2_img_dir + d.bgps2_rms_filen.format(field))
        bgps_img_wcs = _pywcs.WCS(bgps2_img[0].header)
        pix_x, pix_y = _np.floor(bgps_img_wcs.wcs_sky2pix(glon, glat, 0))
        pix_x, pix_y = _np.floor(pix_x)[0], _np.floor(pix_y)[0]
        molcat['map_pk'].ix[cindex] = bgps2_img[0].data[pix_y, pix_x]
        molcat['medmap_pk'].ix[cindex] = bgps2_med[0].data[pix_y, pix_x]
        molcat['map_rms'].ix[cindex] = bgps2_rms[0].data[pix_y, pix_x]
    molcat_pk = molcat[['cnum', 'map_pk', 'medmap_pk', 'map_rms']]
    molcat_pk.to_csv(out_filen + '.csv', index=False, na_rep='-999',
        float_format='%.6f')
    print '-- Catalog written to {}.csv'.format(out_filen)
    return molcat_pk

