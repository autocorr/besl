"""

=========================
Evolutionary Stages Plots
=========================

Plot the BGPS maps and label-masks overplotted with various catalogs for
sign-posts of star-formation and evolutionary stage.

"""
import aplpy
import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib.patheffects as PathEffects
from astropy.io import fits
from besl import catalog
from ..catalog import read_cat
from ..image import get_bgps_img
from ..paths import all_paths as d


def overplot_markers(gc, add_legend=False):
    """
    Overplot markers for different catalogs over an APLpy axis object.
    """
    # ir cats
    ego = read_cat('ego_all')
    msx = read_cat('lumsden13_red_msx')
    msx = msx[msx['Type'].isin(['HII region', 'YSO', 'Young/old star',
                                'HII/YSO'])]
    robit = read_cat('robitaille08_red_spitzer')
    robit_yso = robit[robit['Class'].isin(['cYSO'])]
    robit_agb = robit[robit['Class'].isin(['xAGB', 'sAGB'])]
    hgps = read_cat('higal_70_clean')
    gc.show_markers(hgps['GLON'].values, hgps['GLAT'].values, marker='+',
        color='red', zorder=3, label=r'${\rm HG:70}$')
    gc.show_markers(robit_yso['_Glon'].values, robit_yso['_Glat'].values,
        marker='x', edgecolor='red', facecolor='red', zorder=3,
        label=r'${\rm R08:YSO}$')
    gc.show_markers(robit_agb['_Glon'].values, robit_agb['_Glat'].values,
        marker='x', edgecolor='0.35', facecolor='0.35', zorder=3,
        label=r'${\rm R08:AGB}$')
    gc.show_markers(msx['glon'].values, msx['glat'].values, marker='p',
        edgecolor='red', facecolor='none', zorder=3,
        label=r'${\rm RMS:YSO}$')
    gc.show_markers(ego['_Glon'].values, ego['_Glat'].values, marker='s',
        edgecolor='orange', facecolor='none', zorder=3, label=r'${\rm EGO}$')
    # h2o cats
    h2o_gbt = catalog.read_cat('gbt_h2o')
    h2o_arc = catalog.read_cat('valdettaro01_arcetri')
    h2o_hops = catalog.read_cat('walsh11_hops_h2o')
    gc.show_markers(h2o_gbt['h2o_glon'].values, h2o_gbt['h2o_glat'].values,
        marker='o', edgecolor='green', facecolor='none', zorder=3,
        label=r'${\rm GBT \ H_2O:N}$')
    gc.show_markers(h2o_gbt[h2o_gbt.h2o_f == 1]['h2o_glon'].values,
        h2o_gbt[h2o_gbt.h2o_f == 1]['h2o_glat'].values, marker='o',
        edgecolor='green', facecolor='green', zorder=3,
        label=r'${\rm GBT \ H_2O:Y}$')
    gc.show_markers(h2o_arc[h2o_arc.h2o_f == 1]['_Glon'].values,
        h2o_arc[h2o_arc.h2o_f == 1]['_Glat'].values, marker='^',
        edgecolor='green', facecolor='green', zorder=3,
        label=r'${\rm Lit. \ H_2O}$')
    # TODO fix h2o_rms catalogs
    #h2o_rms = catalog.read_rms_h2o()
    #gc.show_markers(h2o_rms[h2o_rms.h2o_f == 1]['_Glon_x'].values,
    #    h2o_rms[h2o_rms.h2o_f == 1]['_Glat_x'].values, marker='^',
    #    edgecolor='green', facecolor='green', zorder=3)
    gc.show_markers(h2o_hops['lWeight_deg'].values,
        h2o_hops['bWeight_deg'].values, marker='^', edgecolor='green',
        facecolor='green', zorder=3)
    # ch3oh cats
    mmb = read_cat('mmb_all')
    pesta = read_cat('pestalozzi05')
    pandi = read_cat('pandian11')
    gc.show_markers(mmb['glon'].values, mmb['glat'].values, marker='v',
        edgecolor='orange', facecolor='orange', zorder=3,
        label=r'${\rm Lit. CH_3OH}$')
    gc.show_markers(pesta['_Glon'].values, pesta['_Glat'].values, marker='v',
        edgecolor='orange', facecolor='orange', zorder=3)
    gc.show_markers(pandi['glon'].values, pandi['glat'].values, marker='v',
        edgecolor='orange', facecolor='orange', zorder=3)
    # uchii cats
    corn = catalog.read_cornish(exten='hii')
    gc.show_markers(corn['glon'].values, corn['glat'].values, marker='D',
         edgecolor='blue', facecolor='none', zorder=3,
         label=r'${\rm UCHII}$')
    if add_legend:
        _plt.legend(loc='upper center', bbox_to_anchor=(0.125, 0, 0.75, 0.975),
                    bbox_transform=_plt.gcf().transFigure, ncol=5, fontsize=11,
                    scatterpoints=1, markerscale=1.3, frameon=False,
                    mode='expand')
    return gc


def overplot_ellipses(gc):
    """
    Overplot elipses for BGPS sources
    """
    bgps = catalog.read_bgps()
    gc.show_ellipses(bgps['glon_cen'].values, bgps['glat_cen'].values,
        2 * bgps['maj'].values / 3600., 2 * bgps['min'].values / 3600.,
        bgps['pos_ang'].values - 90., edgecolor='grey', facecolor='none',
        alpha=0.75, zorder=2)
    gc.show_markers(bgps['glon_cen'].values, bgps['glat_cen'].values,
        marker='+', edgecolor='grey', facecolor='grey', alpha=0.5, zorder=2)
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
    epsilon = 9e-1 # fudge factor for contouring
    bgps = catalog.read_bgps(exten='all')
    rind = fits.open(d.root_dir + '/BGPS/Images/v2.0.0/' +
                     'v2.0_ds2_{}_13pca_labelmask.fits'.format(field))
    yp, xp = rind[0].data.shape
    X, Y = _np.meshgrid(_np.arange(xp), _np.arange(yp))
    field_cnums = bgps.ix[bgps.field == field, 'field_cnum'].values
    for cnum in field_cnums:
        Z = _np.zeros((yp, xp))
        rind_args = _np.argwhere(rind[0].data == cnum)
        Z[rind_args[:,0] + 1, rind_args[:,1] + 1] = 1
        _plt.contour(X, Y, Z, levels=[epsilon], colors='0.25', linewidths=0.5,
            zorder=1)
    rind.close()
    return gc


def overplot_single_rind(gc, cnum):
    """
    Overplot the label mask contour for a single clump on an APLpy axis object.
    Does not cache between loading fits files, so don't place in loop.

    Parameters
    ----------
    gc : APLpy.FITSFigure
    cnum : number
        BGPS v2.1.0 catalog number

    Returns
    -------
    gc APLpy.FITSFigure
    """
    epsilon = 9e-1  # fudge factor for contouring
    bgps = catalog.read_cat('bgps_v210').set_index('v210cnum')
    field = bgps.loc[cnum, 'field']
    rind = fits.open(d.root_dir + 'BGPS/Images/v2.1.0/v2.1_ds2_' +
                     '{}_13pca_labelmask.fits'.format(field))
    yp, xp = rind[0].data.shape
    X, Y = _np.meshgrid(_np.arange(xp), _np.arange(yp))
    Z = _np.zeros((yp, xp))
    rind_args = _np.argwhere(rind[0].data == cnum)
    Z[rind_args[:,0] + 1, rind_args[:,1] + 1] = 1
    ax = _plt.gca()
    ax.contour(X, Y, Z, levels=[epsilon], colors='k', linewidths=1.0,
               linestyles='solid', zorder=1)
    ax.contour(X, Y, Z, levels=[epsilon], colors='w', linewidths=1.0,
               linestyles='dashed', zorder=2)
    rind.close()
    return gc


def overplot_starless(gc, field):
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
    epsilon = 9e-1 # fudge factor for contouring
    bgps = read_cat('bgps_v210_evo')
    sl_color='#66FF33'
    dets = bgps.query('sf_f != 1 & hg70_eye_f not in [1,2,4]')
    rind = get_bgps_img(field, exten='labelmask', v=210)
    yp, xp = rind[0].data.shape
    X, Y = _np.meshgrid(_np.arange(xp), _np.arange(yp))
    dt_cnums = dets.ix[dets.field == field, 'v210cnum'].values
    for cnum in dt_cnums:
        Z = _np.zeros((yp, xp))
        rind_args = _np.argwhere(rind[0].data == cnum)
        Z[rind_args[:,0] + 1, rind_args[:,1] + 1] = 1
        _plt.contour(X, Y, Z, levels=[epsilon], colors=sl_color,
                     linewidths=2, zorder=1)
    rind.close()
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
    bgps = read_cat('bgps_v210')
    bgps['glon_pix'], bgps['glat_pix'] = gc.world2pixel(bgps['glon_cen'].values,
        bgps['glat_cen'].values)
    bgps['cnum_str'] = bgps['v210cnum'].apply(int).apply(str)
    # get BGPS field cnum's
    bgps_select = bgps.ix[bgps.field == field].reset_index()
    cnums = bgps.ix[bgps.field == field, 'cnum'].values
    # for each source, overplot marker
    for i in xrange(bgps_select.shape[0]):
        txt = _plt.annotate(bgps_select.ix[i, 'cnum_str'],
            xy=(bgps_select.ix[i, 'glon_pix'], bgps_select.ix[i,
            'glat_pix']), xytext=(3.75,-2.5), xycoords='data',
            textcoords='offset points', fontsize='7', weight='book',
            color='black', alpha=0.75, zorder=1)
        txt.set_path_effects([PathEffects.withStroke(linewidth=2,
            foreground='w', alpha=0.75)])
    # TODO overplot green cnum labels for starless clumps
    return gc


def create_stages_plot(lon, lat, dlon, dlat, out_filen='sign_posts'):
    """
    Create two panel figure of BGPS label mask and BGPS image in a specified
    region.

    Parameters
    ----------
    lon : number
        Longitude center coordinate
    lat : number
        Latitude center coordinate
    dlon : number
        Longitude width
    dlat : number
        Latitude height
    out_filen : string, default 'sign_posts'

    Returns
    -------
    gc APLpy.FITSFigure
    """
    ### check if figure is in BGPS
    field_check = catalog.select_bgps_field(lon, lat, coord_type='gal',
            bool_out=True)
    if not field_check:
        raise Exception('Not in BGPS maps.')
    field = catalog.select_bgps_field(lon, lat, coord_type='gal')
    fig = _plt.figure(figsize=(10, 5))
    ### Label mask
    src_img = get_bgps_img(field, exten='labelmask')
    gc1 = aplpy.FITSFigure(src_img,
        figure=fig, subplot=(1,2,1), convention='calabretta')
    gc1.show_colorscale(cmap='gist_gray', vmin=-0.15, vmax=0.05,
        stretch='linear')
    gc1.recenter(lon, lat, width=dlon, height=dlat)
    gc1 = overplot_markers(gc1, add_legend=True)
    gc1 = overplot_rind(gc1, field)
    gc1 = overplot_starless(gc1, field)
    #gc1 = overplot_cnums(gc1, field)
    ### BGPS Image
    flx_img = get_bgps_img(field, exten='map20')
    gc2 = aplpy.FITSFigure(flx_img,
        figure=fig, subplot=(1,2,2), convention='calabretta')
    gc2.show_colorscale(cmap='gist_gray_r', vmin=-0.20, vmax=1.0,
        stretch='linear')
    gc2.recenter(lon, lat, width=dlon, height=dlat)
    gc2.tick_labels.hide_y()
    gc2.axis_labels.hide_y()
    gc2 = overplot_markers(gc2)
    ### Colorbar
    #gc2.add_colorbar(location='right')
    gc_list = [gc1, gc2]
    for gc in gc_list:
        ### Overplot
        ### Grid and ticks
        gc.ticks.set_color('black')
        gc.set_tick_xspacing(0.1)
        gc.set_tick_yspacing(0.1)
        gc.set_tick_labels_xformat('dd.d')
        gc.set_tick_labels_yformat('dd.d')
        gc.add_beam(facecolor='yellow', linewidth=2, hatch='//', major=0.00917,
                    minor=0.00917, angle=0)
    ### Save
    gc2.save(out_filen + '.pdf', dpi=300)
    gc2.save(out_filen + '.eps', dpi=300)
    gc2.save(out_filen + '.png', dpi=300)
    print '-- Printed to file {}'.format(out_filen)
    return gc2


def test_plot():
    #create_stages_plot(23.377173, -0.23825809, 0.5, 0.5)
    #create_stages_plot(32.0, 0, 0.5, 0.5)
    create_stages_plot(12.80, -0.2, 0.5, 0.5)


