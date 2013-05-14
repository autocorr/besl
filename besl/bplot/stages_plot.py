"""

=========================
Evolutionary Stages Plots
=========================

Plot the BGPS maps and label-masks overplotted with various catalogs for
sign-posts of star-formation and evolutionary stage.

"""

import aplpy
import pyfits
import matplotlib.pyplot as _plt
from besl import catalog

def overplot_markers(gc):
    """
    Overplot markers for different catalogs over an APLpy axis object.
    """
    # ir cats
    ego = catalog.read_ego()
    msx = catalog.read_msx()
    robit = catalog.read_robitaille()
    gc.show_markers(ego['_Glon'].values, ego['_Glat'].values, marker='s',
        edgecolor='orange', facecolor='none')
    gc.show_markers(msx['glon'].values, msx['glat'].values, marker='+', color='red')
    gc.show_markers(robit['_Glon'].values, robit['_Glat'].values, marker='x', color='red')
    # h2o cats
    h2o_gbt = catalog.read_gbt_h2o()
    h2o_arc = catalog.read_arcetri_valdettaro()
    h2o_rms = catalog.read_rms_h2o()
    h2o_hops = catalog.read_hops()
    gc.show_markers(h2o_gbt['h2o_glon'].values, h2o_gbt['h2o_glat'].values,
        marker='o', edgecolor='green', facecolor='none')
    gc.show_markers(h2o_gbt[h2o_gbt.h2o_f == 1]['h2o_glon'].values,
        h2o_gbt[h2o_gbt.h2o_f == 1]['h2o_glat'].values, marker='o',
        edgecolor='green', facecolor='green')
    gc.show_markers(h2o_arc[h2o_arc.h2o_f == 1]['_Glon'].values,
        h2o_arc[h2o_arc.h2o_f == 1]['_Glat'].values, marker='^',
        edgecolor='green', facecolor='green')
    gc.show_markers(h2o_rms[h2o_rms.h2o_f == 1]['_Glon_x'].values,
        h2o_rms[h2o_rms.h2o_f == 1]['_Glat_x'].values, marker='^',
        edgecolor='green', facecolor='green')
    gc.show_markers(h2o_hops['lWeight_deg'].values,
        h2o_hops['bWeight_deg'].values, marker='^', edgecolor='green',
        facecolor='green')
    # ch3oh cats
    mmb = catalog.read_mmb()
    gc.show_markers(mmb['glon'].values, mmb['glat'].values, marker='v',
        edgecolor='green', facecolor='green')
    # uchii cats
    corn = catalog.read_cornish(exten='hii')
    gc.show_markers(corn['glon'].values, corn['glat'].values, marker='D',
         edgecolor='blue', facecolor='none')
    return gc

def overplot_ellipses(gc):
    """
    Overplot elipses for BGPS sources
    """
    bgps = catalog.read_bgps()
    gc.show_ellipses(bgps['glon_cen'].values, bgps['glat_cen'].values,
        2 * bgps['maj'].values / 3600., 2 * bgps['min'].values / 3600.,
        bgps['pos_ang'].values - 90., edgecolor='grey', facecolor='none',
        alpha=0.75)
    gc.show_markers(bgps['glon_cen'].values, bgps['glat_cen'].values,
        marker='+', edgecolor='grey', facecolor='grey', alpha=0.5)
    return gc

def overplot_rind(gc, field):
    """
    Overplot label mask contour on an APLpy axis object.

    Parameters
    ----------
    gc : APLpy.FITSFigure
    field : string

    Returns
    -------
    gc APLpy.FITSFigure
    """
    bgps = catalog.read_bgps(exten='all')
    bgps['glon_pix'], bgps['glat_pix'] = gc.world2pixel(bgps['glon'].values,
        bgps['glat'].values)
    rind = pyfits.open('/mnt/eld_data/BGPS/Images/v2.0.0/v2.0_ds2_{}_13pca_labelmask.fits'.format(field))
    # get BGPS field cnum's
    cnums = bgps.ix[bgps.field == field, 'cnum'].values
    for cnum in cnums:
        # for each source, overplot contour
        cnum_rind = rind
        gc.show_contour
    # TODO overplot green contours for starless clumps
    return gc

def overplot_cnums(gc, field):
    """
    Overplot catalog number labels on an APLpy axis object.

    Parameters
    ----------
    gc : APLpy.FITSFigure
    field : string

    Returns
    -------
    gc APLpy.FITSFigure
    """
    bgps = besl.catalog.read_bgps(exten='all')
    bgps['glon_pix'], bgps['glat_pix'] = gc.world2pixel(bgps['glon'].values,
        bgps['glat'].values)
    bgps['cnum_str'] = bgps['cnum'].apply(str)
    rind = pyfits.open('/mnt/eld_data/BGPS/Images/v2.0.0/v2.0_ds2_{}_13pca_labelmask.fits'.format(field))
    # get BGPS field cnum's
    bgps_select = bgps.ix[bgps.field == field].reset_index()
    cnums = bgps.ix[bgps.field == field, 'cnum'].values
    # for each source, overplot marker
    for i in xrange(bgps_select.shape[0]):
        txt = _plt.annotate(bgps.ix[i, 'cnum_str'], xy=(bgps.ix[i, 'glon_pix'],
            bgps.ix[i, 'glat_pix']), xytext=(3.75,-2.5), xycoords='data',
            textcoords='offset points', fontsize='6', weight='book',
            color='grey')
        txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w')])
    # TODO overplot green cnum labels for starless clumps
    return gc

def create_stages_plot(lon, lat, dlon, dlat, out_filen='sign_posts'):
    # check if figure is in BGPS
    field_check = catalog.select_bgps_field(lon, lat, coord_type='gal',
            bool_out=True)
    if not field_check:
        raise Exception('Not in BGPS maps.')
    field = catalog.select_bgps_field(lon, lat, coord_type='gal')
    fig = _plt.figure(figsize=(10, 5))
    # Label mask
    gc1 = aplpy.FITSFigure(
        '/mnt/eld_data/BGPS/Images/v2.0.0/v2.0_ds2_{}_13pca_labelmask.fits'.format(field),
        figure=fig, subplot=(1,2,1), convention='calabretta')
    gc1.show_colorscale(cmap='gist_gray', vmin=-0.15, vmax=0.05,
        stretch='linear')
    gc1.recenter(lon, lat, width=dlon, height=dlat)
    gc1 = overplot_rind(gc1, field)
    gc1 = overplot_cnums(gc1, field)
    # BGPS Image
    gc2 = aplpy.FITSFigure(
        '/mnt/eld_data/BGPS/Images/v2.0.0/v2.0_ds2_{}_13pca_map20.fits'.format(field),
        figure=fig, subplot=(1,2,2), convention='calabretta')
    gc2.show_colorscale(cmap='gist_gray_r', vmin=-0.20, vmax=1.0,
        stretch='linear')
    gc2.recenter(lon, lat, width=dlon, height=dlat)
    gc2.hide_ytick_labels()
    gc2.hide_yaxis_label()
    # Colorbar
    #gc2.add_colorbar(location='right')
    gc_list = [gc1, gc2]
    for gc in gc_list:
        # Overplot
        gc = overplot_markers(gc)
        #gc = overplot_ellipses(gc)
        # Grid and ticks
        gc.ticks.set_color('black')
        gc.set_tick_xspacing(0.1)
        gc.set_tick_yspacing(0.1)
        gc.set_tick_labels_xformat('dd.d')
        gc.set_tick_labels_yformat('dd.d')
        gc.add_beam(facecolor='yellow', linewidth=2, hatch='//', major=0.00917,
            minor=0.00917, angle=0)
    gc.save(out_filen + '.eps', dpi=300)
    gc.save(out_filen + '.png', dpi=300)
    print '-- Printed to file {}'.format(out_filen)
    return gc

def test_plot():
    create_stages_plot(23.377173, -0.23825809, 0.5, 0.5)

