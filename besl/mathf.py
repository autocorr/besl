"""
==============
Math Functions
==============

Mathematical functions.

"""

import numpy as _np
from scipy.special import ellipe


def ang_diff(a, b, angle_type='deg'):
    """
    Take the difference of two angles between 0 - 360 degrees.

    Parameters
    ----------
    a, b : number
        Numbers representing two angles.
    angle_type : str, default 'deg'
        Type of angle, valid types: 'deg' or 'rad'.

    Returns
    -------
    number : number
    """
    if angle_type == 'deg':
        return 180 - _np.abs(_np.abs(a - b) - 180)
    elif angle_type == 'rad':
        return _np.pi / 2. - _np.abs(_np.abs(a - b) - _np.pi / 2.)
    else:
        raise ValueError('Invalid angle type: {0}.'.format(angle_type))


def bin_factor(x):
    """
    Factor a decimal integer into it's power of 2 components.

    Parameters
    ----------
    x : int
        Must be a positive integer

    Returns
    -------
    factors : numpy.array
        Array of non-zero binary factors
    """
    if (not isinstance(x, int)) | (x < 0):
        raise ValueError('Must be a positive integer: {0}.'.format(x))
    b = _np.array(list(bin(x)[2:]))[::-1].astype(int)
    factors = (1 << _np.arange(len(b)))[b.astype(bool)]
    return factors


def weighted_mean(a, a_err):
    """
    Compute the weighted mean.

    Parameters
    ----------
    a : array_like
        Values.
    a_err : array_like
        Uncertainties of values in a. Use masked array if
        array contains zero valued elements.
    """
    return _np.sum(a / a_err**2) / _np.sum(1 / a_err**2)


def med_abs_dev(a):
    """
    Compute the median absolute deviation.

    Parameters
    ----------
    a : array_like
        Values.
    """
    return _np.median(_np.abs(a - _np.median(a)))


def sigfigs(a, figs=1):
    """
    Rounds values to fixed number of significant digits.

    Parameters
    ----------
    a : array_like
        Values.
    figs : number
        Number of significant figures
    """
    figs = _np.int(figs)
    a = _np.array([_np.round(x, -_np.int(_np.floor(_np.log10(x))) + figs - 1) for x
        in a])
    return a


def arcsinhspace(start, stop, num=50):
    """
    Return numbers spaced evenly on an arcsinh scale.

    Parameters
    ----------
    start : number
        start value
    stop : number
        stop/end value, inclusive
    num : number
        number of intervales between start and stop
    """
    return _np.sinh(_np.linspace(_np.arcsinh(start), _np.arcsinh(stop), num))


def ellipse_area(major, minor, semi=True):
    """
    Calculate the area of an ellipse given the axe lengths.

    Parameters
    ----------
    major : number
        Major axis length
    minor : number
        Minor axis length
    semi : Bool, default True
        Use semi- (True) or full-axes (False) lengths

    Returns
    -------
    area : number
        Area of ellipse
    """
    if semi:
        return _np.pi * major * minor
    elif not semi:
        return _np.pi * major * minor / 4.


def ellipse_circumf(major, minor, semi=True):
    """
    Calculate the circumference of an ellipse given the axe lengths.

    Parameters
    ----------
    major : number
        Major axis length
    minor : number
        Minor axis length
    semi : Bool, default True
        Use semi- (True) or full-axes (False) lengths

    Returns
    -------
    area : number
        Area of ellipse
    """
    eccen2 = (major**2 - minor**2) / major**2
    if semi:
        return 4. * major * ellipe(eccen2)
    elif not semi:
        return 2. * major * ellipe(eccen2)


def pd_weighted_mean(df, val_cols, err_cols, out_cols=('val', 'err')):
    """
    Weighted mean of a DataFrame for multiple columns in both values and
    associated uncertainties. Inverse uncertainties are used as weights. The
    full DataFrame is returned with columns added for weighted values and
    weighted errors.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing value and error columns
    val_cols : list-like
        Column names for values
    err_cols : list-like
        Column names for associated errors in the same order as `val_cols`.
    out_cols : tuple, default ('val', 'err')
        Default names for weighted columns

    Returns
    -------
    df : pd.DataFrame
        Full DataFrame with weighted column names added
    """
    assert len(val_cols) > 1
    assert len(err_cols) > 1
    errs_to_vals = {e: v for v, e in zip(val_cols, err_cols)}
    # Sum of weights
    sum_weights = (1. / df[err_cols]).sum(axis=1)
    # Weighted values
    # note: `.div` automatically divides columns of the same name
    #       and axis 1 is by row
    weighted_vals = df[val_cols].div(
                    df[err_cols].rename(columns=errs_to_vals)).sum(axis=1)
    df[out_cols[0]] = weighted_vals.div(sum_weights)
    # Weighted errors by adding in quadrature
    nerrs = df[err_cols].count(axis=1).astype(float)
    nerrs[nerrs == 0] = _np.nan
    df[out_cols[1]] = (_np.sqrt(df[err_cols].pow(2).sum(axis=1))).div(nerrs)
    return df


