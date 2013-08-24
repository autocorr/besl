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
    if angle_type='deg':
        return 180 - _np.abs(_np.abs(a - b) - 180)
    elif angle_type='rad':
        return _np.pi / 2. - _np.abs(_np.abs(a - b) - _np.pi / 2.)
    else:
        raise ValueError('Invalid angle type: {0}.'.format(angle_type))

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


