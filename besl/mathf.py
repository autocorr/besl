"""
Mathematical functions.
"""

import numpy as _np

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

