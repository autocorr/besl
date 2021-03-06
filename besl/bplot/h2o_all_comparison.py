"""
=======================
H2O All Comparison Plot
=======================

Plotting routines for comparing BGPS GBT, HOPS, RMS, and Arectri H2O maser
surveys and catalogs.

"""

import numpy as _np
import matplotlib.pyplot as _plt

def h2o_sample_statistics():
    """
    Print medians and median absolute deviations for H2O maser observations.
    Includes BGPS GBT, RMS H2O, Arcetri (Valdettaro), and HOPS
    """
    from besl.catalog import read_gbt_h2o, read_arcetri_cesaroni, \
        read_arcetri_valdettaro, read_hops, read_rms_h2o, match_h2o_maser
    from besl.mathf import med_abs_dev
    # read in catalogs
    gbt_h2o = read_gbt_h2o()
    rms_h2o = read_rms_h2o()
    arc_val = read_arcetri_valdettaro()
    hops = read_hops()
    gbt_h2o_det = gbt_h2o[gbt_h2o.h2o_f == 1]
    arc_val_det = arc_val[arc_val.h2o_f == 1]
    # print stats
    print '-- median +/- median absolute deviation'
    print '-----------------------------------------------------'
    print '-- BGPS GBT'
    print '-- detections: {}'.format(gbt_h2o_det.shape[0])
    print '-- peak intensity [Jy]: {} +/- {}'.format(
        _np.median(gbt_h2o_det.h2o_tpk_jy),
        med_abs_dev(gbt_h2o_det.h2o_tpk_jy))
    print '-- peak intensity uncertainty [Jy]: {} +/- {}'.format(
        _np.median(gbt_h2o_det.h2o_tpk_err_jy),
        med_abs_dev(gbt_h2o_det.h2o_tpk_err_jy))
    print '-- integrated intensity [Jy km/s]: {} +/- {}'.format(
        _np.median(gbt_h2o_det.h2o_int_jy),
        med_abs_dev(gbt_h2o_det.h2o_int_jy))
    print '-- integrated intensity uncertainty [Jy km/s]: {} +/- {}'.format(
        _np.median(gbt_h2o_det.h2o_int_err_jy),
        med_abs_dev(gbt_h2o_det.h2o_int_err_jy))
    print '-- velocity spread [km/s]: {} +/- {}'.format(
        _np.median(gbt_h2o_det.h2o_vsp),
        med_abs_dev(gbt_h2o_det.h2o_vsp))
    print '-----------------------------------------------------'
    print '-- RMS H2O'
    print '-- detections: {}'.format(rms_h2o.shape[0])
    print '-- peak intensity [Jy]: {} +/- {}'.format(
        _np.median(rms_h2o.Fpeak),
        med_abs_dev(rms_h2o.Fpeak))
    print '-- velocity spread [km/s]: {} +/- {}'.format(
        _np.median(rms_h2o.vsp),
        med_abs_dev(rms_h2o.vsp))
    print '-----------------------------------------------------'
    print '-- Arcetri (Valdettaro)'
    print '-- detections: {}: '.format(arc_val_det.shape[0])
    print '-- peak intensity [Jy]: {} +/- {}'.format(
        _np.median(arc_val_det.Fp),
        med_abs_dev(arc_val_det.Fp))
    print '-- peak intensity uncertainty [Jy]: {} +/- {}'.format(
        _np.median(arc_val_det.Sig),
        med_abs_dev(arc_val_det.Sig))
    print '-- integrated intensity [Jy]: {} +/- {}'.format(
        _np.median(arc_val_det.Stot),
        med_abs_dev(arc_val_det.Stot))
    print '-- velocity spread [km/s]: {} +/- {}'.format(
        _np.median(arc_val_det.vsp),
        med_abs_dev(arc_val_det.vsp))
    print '-----------------------------------------------------'
    print '-- HOPS'
    print '-- detections: '.format(hops.shape[0])
    print '-- peak intensity [Jy]: {} +/- {}'.format(
        _np.median(hops.T_peak_Jy),
        med_abs_dev(hops.T_peak_Jy))
    print '-- peak intensity uncertainty [Jy]: {} +/- {}'.format(
        _np.median(hops.RMS_Jy),
        med_abs_dev(hops.RMS_Jy))
    print '-- integrated intensity [Jy km/s]: {} +/- {}'.format(
        _np.median(hops.Tbdv_Jykms),
        med_abs_dev(hops.Tbdv_Jykms))
    print '-- integrated intensity uncertainty [Jy km/s]: {} +/- {}'.format(
        _np.median(hops.dTbdv_Jykms),
        med_abs_dev(hops.dTbdv_Jykms))
    print '-- velocity spread [km/s]: {} +/- {}'.format(
        _np.median(hops.vsp),
        med_abs_dev(hops.vsp))
    print '-----------------------------------------------------'
    return

def all_maser_peak_intensity_comparison():
    """
    Plot histograms of H2O maser observations peak intensities and integrated
    intensities.  Includes BGPS GBT, RMS H2O, Arcetri (Valdettaro), and HOPS

    Returns
    -------
    fig : matplotlib figure object
    ax1 : matplotlib axes object
        Intensity histogram sorted by catalog
    ax2 : matplotlib axes object
        Integrated intensity histogram sorted by catalog
    """
    import os
    from besl.catalog import read_gbt_h2o, read_arcetri_valdettaro, \
        read_hops, read_rms_h2o, match_h2o_maser
    # read in catalogs
    gbt_h2o = read_gbt_h2o()
    rms_h2o = read_rms_h2o()
    arc_val = read_arcetri_valdettaro()
    hops = read_hops()
    gbt_h2o_det = gbt_h2o[gbt_h2o.h2o_f == 1]
    arc_val_det = arc_val[arc_val.h2o_f == 1]
    # plot settings
    lbins_pk = _np.logspace(_np.log10(0.1), _np.log10(1e4), 40)
    lbins_int = _np.logspace(_np.log10(0.05), _np.log10(1e6), 40)
    hist_kwargs_gbt = {'label': 'BGPS GBT', 'color': 'green',
            'histtype': 'step', 'linestyle': 'solid', 'linewidth': 2, 'log':
            True, 'hatch': '.'}
    hist_kwargs_rms = {'label': 'RMS', 'color': 'red',
            'histtype': 'step', 'hatch': '/', 'linestyle': 'dashed',
            'linewidth': 2, 'log': True}
    hist_kwargs_arc = {'label': 'Arcetri', 'color': 'lightslategray',
            'histtype': 'step', 'hatch': '\\', 'linestyle': 'solid',
            'linewidth': 2, 'log': True}
    hist_kwargs_hop = {'label': 'HOPS', 'color': 'blue',
            'histtype': 'step', 'linestyle': 'dashed', 'linewidth': 2, 'log':
            True}
    fig = _plt.figure(figsize=(8, 6))
    # begin plot
    ax1 = fig.add_subplot(211)
    ax1.hist(gbt_h2o_det.h2o_tpk_jy, bins=lbins_pk, **hist_kwargs_gbt)
    ax1.hist(rms_h2o.Fpeak, bins=lbins_pk, **hist_kwargs_rms)
    ax1.hist(arc_val.Fp, bins=lbins_pk, **hist_kwargs_arc)
    ax1.hist(hops.T_peak_Jy, bins=lbins_pk, **hist_kwargs_hop)
    ax1.set_xlabel(r'$I \ \ [{\rm Jy}]$')
    ax1.set_ylabel(r'$N$')
    ax1.set_xlim([8e-2, 1e4])
    ax1.set_ylim([8e-1, 3e2])
    ax1.set_xscale('log')
    ax1.legend(loc=0, prop={'size':12})
    ax2 = fig.add_subplot(212)
    ax2.hist(gbt_h2o_det.h2o_int_jy, bins=lbins_int, **hist_kwargs_gbt)
    ax2.hist(arc_val.Stot, bins=lbins_int, **hist_kwargs_arc)
    ax2.hist(hops.Tbdv_Jykms, bins=lbins_int, **hist_kwargs_hop)
    ax2.set_xlabel(r'$\int I \ {\rm d}v \ \ [{\rm Jy \ km \ s^{-1}}]$')
    ax2.set_ylabel(r'$N$')
    ax2.set_xlim([5e-2, 1e6])
    ax2.set_ylim([8e-1, 2e2])
    ax2.set_xscale('log')
    ax2.legend(loc=0, prop={'size':12})
    _plt.show()
    _plt.savefig(os.getcwd() + '/h2o_all_peak_int_hists.png', format='png',
        dpi=200)
    return [fig, ax1, ax2]

def rms_classification_comparison():
    """
    Plot histograms of H2O maser observations in the RMS catalog by
    classification.

    Returns
    -------
    fig : matplotlib figure object
    ax1 : matplotlib axes object
        Velocity spread histogram sorted by catalog
    """
    import os
    from besl.catalog import read_gbt_h2o, read_arcetri_valdettaro, \
        read_hops, read_rms_h2o, match_h2o_maser
    # read in catalogs
    rms_h2o = read_rms_h2o()
    # plot settings
    lbins_pk = _np.logspace(_np.log10(0.1), _np.log10(1e3), 20)
    hist_kwargs_yso = {'label': 'YSO', 'color': 'green',
            'histtype': 'step', 'linestyle': 'solid', 'linewidth': 2, 'log':
            True}
    hist_kwargs_hii = {'label': 'HII region', 'color': 'blue',
            'histtype': 'step', 'linestyle': 'dashed',
            'linewidth': 2, 'log': True}
    hist_kwargs_hiiyso = {'label': 'HII/YSO', 'color': 'cyan',
            'histtype': 'step', 'linestyle': 'solid',
            'linewidth': 2, 'log': True}
    hist_kwargs_star = {'label': 'Young/old star', 'color': 'red',
            'histtype': 'step', 'linestyle': 'dashed', 'linewidth': 2, 'log':
            True}
    hist_kwargs_unk = {'label': 'Unknown', 'color': 'gray',
            'histtype': 'step', 'linestyle': 'dashed', 'linewidth': 2, 'log':
            True}
    fig = _plt.figure(figsize=(8, 3))
    # begin plot
    ax1 = fig.add_subplot(111)
    ax1.hist(rms_h2o[rms_h2o.Type == 'YSO'].Fpeak, bins=lbins_pk,
        **hist_kwargs_yso)
    #ax1.hist(rms_h2o[rms_h2o.Type == 'HII region'].Fpeak, bins=lbins_pk,
    #    **hist_kwargs_hii)
    #ax1.hist(rms_h2o[rms_h2o.Type == 'HII/YSO'].Fpeak, bins=lbins_pk,
    #    **hist_kwargs_hiiyso)
    ax1.hist(rms_h2o[rms_h2o.Type == 'Young/old star'].Fpeak, bins=lbins_pk,
        **hist_kwargs_star)
    #ax1.hist(rms_h2o[rms_h2o.Type.isnull()].Fpeak, bins=lbins_pk,
    #    **hist_kwargs_unk)
    ax1.set_xlabel(r'$I \ \ [{\rm Jy}]$')
    ax1.set_ylabel(r'$N$')
    ax1.set_xlim([8e-2, 1e4])
    ax1.set_ylim([8e-1, 2e1])
    ax1.set_xscale('log')
    ax1.legend(loc=0, prop={'size':12})
    _plt.show()
    _plt.savefig(os.getcwd() + '/rms_classification_comparison.pdf',
        format='pdf', dpi=200)
    return [fig, ax1]

def all_maser_width_comparison():
    """
    Plot histograms of H2O maser observations.  Includes BGPS GBT, RMS H2O,
    Arcetri (Valdettaro), and HOPS

    Returns
    -------
    fig : matplotlib figure object
    ax1 : matplotlib axes object
        Intensity histogram
    """
    import os
    from besl.catalog import read_gbt_h2o, read_arcetri_valdettaro, \
        read_hops, read_rms_h2o, match_h2o_maser
    # read in catalogs
    gbt_h2o = read_gbt_h2o()
    rms_h2o = read_rms_h2o()
    arc_val = read_arcetri_valdettaro()
    hops = read_hops()
    gbt_h2o_det = gbt_h2o[gbt_h2o.h2o_f == 1]
    arc_val_det = arc_val[arc_val.h2o_f == 1]
    # plot settings
    lbins_vsp = _np.logspace(_np.log10(0.1), _np.log10(1e3), 40)
    hist_kwargs_gbt = {'bins': lbins_vsp, 'label': 'BGPS GBT', 'color': 'green',
            'histtype': 'step', 'linestyle': 'solid', 'linewidth': 2, 'log':
            True, 'hatch': '.'}
    hist_kwargs_rms = {'bins': lbins_vsp, 'label': 'RMS', 'color': 'red',
            'histtype': 'step', 'hatch': '/', 'linestyle': 'dashed',
            'linewidth': 2, 'log': True}
    hist_kwargs_arc = {'bins': lbins_vsp, 'label': 'Arcetri', 'color': 'lightslategray',
            'histtype': 'step', 'hatch': '\\', 'linestyle': 'solid',
            'linewidth': 2, 'log': True}
    hist_kwargs_hop = {'bins': lbins_vsp, 'label': 'HOPS', 'color': 'blue',
            'histtype': 'step', 'linestyle': 'dashed', 'linewidth': 2, 'log':
            True}
    fig = _plt.figure(figsize=(8, 6))
    # begin plot
    ax1 = fig.add_subplot(111)
    ax1.hist(gbt_h2o_det.h2o_vsp, **hist_kwargs_gbt)
    ax1.hist(rms_h2o.vsp, **hist_kwargs_rms)
    ax1.hist(arc_val.vsp, **hist_kwargs_arc)
    ax1.hist(hops.vsp, **hist_kwargs_hop)
    ax1.set_xlabel(r'$v_{\rm spread} \ \ [{\rm km \ s^{-1}}]$')
    ax1.set_ylabel(r'$N$')
    ax1.set_xlim([1e-1, 1e3])
    ax1.set_ylim([8e-1, 3e2])
    ax1.set_xscale('log')
    ax1.legend(loc=9, prop={'size':12})
    _plt.show()
    _plt.savefig(os.getcwd() + '/h2o_all_width_hists.png', format='png',
        dpi=200)
    return [fig, ax1]

