"""
==================
BGPS Finder Charts
==================

Create BGPS finder charts from the v2 maps.

"""

import ipdb as pdb
import aplpy
import pyfits
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.patheffects as PathEffects
from besl import catalog

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
dets = molcat[molcat.mol_vlsr_f.isin([1,3])]
mdet = molcat[molcat.mol_vlsr_f == 2]
ndet = molcat[molcat.mol_vlsr_f == 0]
img_list = bounds[bounds.glon_min > 7].field.values
# Calculate round coordinate intervals
bounds['glon_round_sep'] = (bounds['glon_max'] - bounds['glon_min']) / \
    _np.round(bounds['glon_max'] - bounds['glon_min'])
bounds['glat_round_sep'] = (bounds['glat_max'] - bounds['glat_min']) / \
    _np.round(bounds['glat_max'] - bounds['glat_min'])
# Sort by Galactic longitude
bounds = bounds.sort('glon_min', ascending=True)
# String labels
dets['mol_vlsr_str'] = dets.mol_vlsr_f.apply(str)
molcat['cnum_str'] = molcat.cnum.apply(str)

def create_field_tile(img_fits, field='temp', show_cnum=False):
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
    ### Plot
    fig = _plt.figure(figsize=(15,7))
    gc = aplpy.FITSFigure(img_fits, figure=fig, convention='calabretta')
    gc.show_colorscale(cmap='gist_gray_r', vmin=-0.20, vmax=0.65, stretch='linear')
    ### Cnum labels
    if show_cnum:
        molcat['glon_pix'], molcat['glat_pix'] = gc.world2pixel(
            molcat.hht_glon.values, molcat.hht_glat.values)
        for i in xrange(molcat.shape[0]):
            txt = _plt.annotate(molcat.cnum_str.values[i],
                xy=(molcat.glon_pix.values[i], molcat.glat_pix.values[i]),
                xytext=(3.75,-2.5), xycoords='data', textcoords='offset points',
                fontsize='7', weight='book')
            txt.set_path_effects([PathEffects.withStroke(linewidth=2,
                foreground='w')])
    ### Velocity labels
    dets['glon_pix'], dets['glat_pix'] = gc.world2pixel(dets.hht_glon.values,
        dets.hht_glat.values)
    for i in xrange(dets.shape[0]):
        txt = _plt.annotate(dets.mol_vlsr_str.values[i],
            xy=(dets.glon_pix.values[i], dets.glat_pix.values[i]),
            xytext=(3.75,3.5), xycoords='data', textcoords='offset points',
            fontsize='8', weight='book', color='green')
        txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
    ### Markers
    gc.show_markers(dets.hht_glon.values, dets.hht_glat.values,
        layer='detects', marker='o', linewidth=2, edgecolor='green')
    gc.show_markers(mdet.hht_glon.values, mdet.hht_glat.values,
        layer='multi_detects', marker='d', linewidth=2, edgecolor='yellow')
    gc.show_markers(ndet.hht_glon.values, ndet.hht_glat.values,
        layer='non_detects', marker='x', linewidth=2, edgecolor='red')
    ### Grid and ticks
    gc.ticks.set_color('black')
    gc.set_tick_xspacing(0.1)
    gc.set_tick_yspacing(0.1)
    gc.set_tick_labels_xformat('dd.d')
    gc.set_tick_labels_yformat('dd.d')
    gc.show_grid()
    ### Beam and colorbar
    gc.add_beam(facecolor='yellow', linewidth=2, hatch='//', major=0.00917,
            minor=0.00917, angle=0)
    # TODO add label to colorbar
    gc.add_colorbar(location='top')#, xlabel=r'${\rm Jy \ beam^{-1}}$')
    print '-- Image Loaded: {0}'.format(field)
    return gc

def create_sub_tile(gc, field, i, lon_pos, lat_pos, dlon, dlat):
    """
    Recenter and zoom on a BGPS field to create a sub-plot finder chart from a
    APLpy.Axes object.

    Paramters
    ---------
    gc : APLpy.Axes
        BGPS field axes object
    field : string
        Field name to tag the output file
    i : number
        Absolute chart index number
    lon_pos, lat_pos : number
        Galactic longitude and latitude in decimal degrees
    dlon, dlat : number
        Chart width and height in decimal degrees
    """
    gc.recenter(lon_pos, lat_pos, width=dlon, height=dlat)
    gc.save('{0}{1:0>3d}_{2}_{3:.2f}{4:+.2f}.eps'.format(d.out_dir, i, field,
        lon_pos, lat_pos), dpi=250)
    gc.save('{0}{1:0>3d}_{2}_{3:.2f}{4:+.2f}.png'.format(d.out_dir, i, field,
        lon_pos, lat_pos), dpi=250)
    print '-- Printing: {0} {1:.2f} {2:+.2f}'.format(field, lon_pos, lat_pos)
    return

def is_in_subtile(lon_cen, lat_cen, dlon, dlat):
    """
    Return whether there are any BGPS sources within a range. Input values in
    decimal degrees from 0 to 360.

    Parameters
    ----------
    lon_cen : number
        Longitude center position
    lat_cen : number
        Latitude center position
    dlon : number
        Width of tile
    dlat : number
        Height of tile

    Returns
    -------
    Boolean
    """
    # TODO fix for sources near 0 - 360 boundary
    num_in_tile = molcat[(molcat.hht_glon > lon_cen - dlon / 2.) &
                         (molcat.hht_glon < lon_cen + dlon / 2.) &
                         (molcat.hht_glat > lat_cen - dlat / 2.) &
                         (molcat.hht_glat < lat_cen + dlat / 2.)].shape[0]
    if num_in_tile > 0:
        return True
    elif num_in_tile == 0:
        return False
    else:
        raise Exception

def create_velo_finder_charts():
    """
    Generate finder charts for the BGPS Molecular Line survey.
    """
    cnum_list = open('chart_source_list.cat', 'w')
    cnum_list.write('cnum,chart_index\n')
    img_list = bounds[(bounds.glon_min > 7.4) &
                      (bounds.glon_max < 197)].field.values
    # Select fields and print tiles
    i = 1
    for field in img_list:
        img_fits = pyfits.open(d.img_dir + d.img_filen.format(field))
        gc = create_field_tile(img_fits, field, show_cnum=False)
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
                if is_in_subtile(lon_pos, lat_pos, dlon, dlat):
                    create_sub_tile(gc, field, i, lon_pos, lat_pos,
                        dlon=dlon, dlat=dlat)
                    subtile_sources = molcat[
                            (molcat.hht_glon > lon_pos - dlon / 2.) &
                            (molcat.hht_glon < lon_pos + dlon / 2.) &
                            (molcat.hht_glat > lat_pos - dlat / 2.) &
                            (molcat.hht_glat < lat_pos + dlat / 2.)]
                    for in_tile_cnum in subtile_sources.cnum:
                        cnum_list.write('{0:0<4d},{1:0>3d}\n'.format(in_tile_cnum,i))
                    i += 1
        gc.close()
        img_fits.close()
    cnum_list.close()
    return

