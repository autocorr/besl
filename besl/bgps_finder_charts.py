"""
==================
BGPS Finder Charts
==================

Create BGPS finder charts from the v2 maps.

"""

import aplpy
import pyfits
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.patheffects as PathEffects
import catalog

class Dirs(object):
    def __init__(self):
        self.img_dir = '/mnt/eld_data/BGPS/Images/v2.0.0/'
        self.img_filen = 'v2.0_ds2_{}_13pca_map20.fits'
        self.out_dir = '/mnt/eld_data/BGPS/molec_finding_charts/out_plots/v2.0.0/'
d = Dirs()

# Read in data
bgps = catalog.read_bgps()
molcat = catalog.read_molcat()
bounds = catalog.read_bgps_bounds()
# HCO+ flag groups
dets = molcat[(molcat.hco_f != 0) & (molcat.hco_f != 2)]
mdet = molcat[molcat.hco_f == 2]
ndet = molcat[molcat.hco_f == 0]
img_list = bounds[bounds.glon_min > 7].field.values
# Calculate round coordinate intervals
bounds['glon_round_sep'] = (bounds['glon_max'] - bounds['glon_min']) / \
    _np.round(bounds['glon_max'] - bounds['glon_min'])
bounds['glat_round_sep'] = (bounds['glat_max'] - bounds['glat_min']) / \
    _np.round(bounds['glat_max'] - bounds['glat_min'])

def create_field_tile(img_fits, field='temp'):
    """
    Create a BGPS finder chart for a field at specified coordinates.

    Parameters
    ----------
    img_fits : pyfits.hdu
        BGPS image fits to plot
    field : string, default 'temp'
        BGPS field name

    Returns
    -------
    gc : APLpy.Axes
    """
    ### plot
    fig = _plt.figure(figsize=(15,7))
    gc = aplpy.FITSFigure(img_fits, figure=fig, convention='calabretta')
    gc.show_colorscale(cmap='gist_gray_r', vmin=-0.20, vmax=0.65, stretch='linear')
    ### velocity labels
    dets['hco_v_str'] = dets.hco_v.apply(str)
    dets['glon_pix'], dets['glat_pix'] = gc.world2pixel(dets.hht_glon.values,
        dets.hht_glat.values)
    for i in xrange(dets.shape[0]):
        txt = _plt.annotate(dets.hco_v_str.values[i],
            xy=(dets.glon_pix.values[i], dets.glat_pix.values[i]),
            xytext=(3,3), xycoords='data', textcoords='offset points',
            fontsize='10', weight='book')
        txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
    ### markers
    gc.show_markers(dets.hht_glon.values, dets.hht_glat.values,
        layer='detects', marker='o', linewidth=2, edgecolor='green')
    gc.show_markers(mdet.hht_glon.values, mdet.hht_glat.values,
        layer='multi_detects', marker='d', linewidth=2, edgecolor='yellow')
    gc.show_markers(ndet.hht_glon.values, ndet.hht_glat.values,
        layer='non_detects', marker='x', linewidth=2, edgecolor='red')
    ### grid and ticks
    gc.ticks.set_color('black')
    gc.set_tick_xspacing(0.1)
    gc.set_tick_yspacing(0.1)
    gc.set_tick_labels_xformat('dd.d')
    gc.set_tick_labels_yformat('dd.d')
    gc.show_grid()
    ### beam and colorbar
    gc.add_beam(facecolor='yellow', linewidth=2, hatch='//', major=0.00917,
            minor=0.00917, angle=0)
    # TODO add label to colorbar
    gc.add_colorbar(location='top')#, xlabel=r'${\rm Jy \ beam^{-1}}$')
    print '-- Image Loaded: {0}'.format(field)
    return gc

def create_sub_tile(gc, field, lon_pos, lat_pos, dlon, dlat):
    """
    Recenter and zoom on a BGPS field to create a sub-plot finder chart from a
    APLpy.Axes object.

    Paramters
    ---------
    gc : APLpy.Axes
        BGPS field axes object
    field : string
        Field name to tag the output file
    lon_pos, lat_pos : number
        Galactic longitude and latitude in decimal degrees
    dlon, dlat : number
        Chart width and height in decimal degrees
    """
    gc.recenter(lon_pos, lat_pos, width=dlon, height=dlat)
    gc.save('{0}{1}_{2:.2f}{3:+.2f}.png'.format(d.out_dir, field, lon_pos,
        lat_pos), dpi=150)
    print '-- Printing: {0} {1:.2f} {2:+.2f}'.format(field, lon_pos, lat_pos)
    return

def create_velo_finder_charts():
    """
    Generate finder charts for the BGPS Molecular Line survey.
    """
    v2_new_fields = ['camob1', 'iras_22172', 'sh235']
    img_list = bounds[(bounds.glon_min > 90) & (bounds.glon_max < 197) &
            _np.logical_not(bounds.field.isin(v2_new_fields))].field.values
    # Select fields and print tiles
    for field in img_list:
        img_fits = pyfits.open(d.img_dir + d.img_filen.format(field))
        gc = create_field_tile(img_fits, field)
        lon_min, lon_max, lat_min, lat_max = bounds.ix[bounds.field ==
            field][['glon_min', 'glon_max', 'glat_min', 'glat_max']].values[0]
        dlon = 0.5 * \
            bounds.ix[bounds.field == field]['glon_round_sep'].values[0]
        dlat = 0.5 * \
            bounds.ix[bounds.field == field]['glat_round_sep'].values[0]
        lon_ints = _np.arange(lon_min + dlon / 2, lon_max + dlon / 2, dlon)
        lat_ints = _np.arange(lat_min + dlat / 2, lat_max + dlat / 2, dlat)
        for lon_pos in lon_ints:
            for lat_pos in lat_ints:
                create_sub_tile(gc, field, lon_pos, lat_pos,
                    dlon=dlon, dlat=dlat)
        gc.close()
        img_fits.close()
    return

