"""
Routines for manipulating coordinates.
"""

import math as _m
import numpy as _np

def dec2sexstr(deci, sfigs=1, hd='h', lead_psign=False):
    """
    Convert decimal degree to a sexagesimal string of the form 'HH:MM:SS.S' by
    default.

    Parameters
    ----------
    sfigs : number
        Number of significant figures after decimal in seconds
    hd : string, ('h', 'd')
        Hour or degree convention
    lead_sign : Boolean
        Whether to include the leading sign +/- in string
    """
    lead_psymb = ''
    dform = {'h': 24 / 360., 'd': 1.}
    (h_frac, h) = _np.modf(deci * dform[hd])
    (m_frac, m) = _np.modf(h_frac * 60)
    s = m_frac * 60
    if (lead_psign is True) & (h >= 0):
        lead_psymb = '+'
    coord_string = "{0}{1:0>2d}:{2:0>2d}:{3:0>2d}.{4:0>{sfigs}d}".format(
        lead_psymb, int(h), _np.abs(int(m)), _np.abs(int(s)),
        _np.abs(int(_np.mod(s,1) * 10**sfigs)), sfigs=sfigs)
    return coord_string

def sexstr2dec(sexstr, sep=':', hd='h'):
    """
    Convert a sexagesimal string of delimited by a seperator character, eg
    "+HH:MM:SS.S" with ":", into a decimal float.

    Parameters
    ----------
    sep : string
        Seperator character between hours, minutes, and seconds
    hd : string, ('h', 'd')
        Hour or degree convention
    """
    dform = {'h': 24. / 360., 'd': 1.}
    h, m, s = [float(i) for i in sexstr.split(sep)]
    deci = _np.sign(h) * (_np.abs(h) + m / 60. + s / 3600.) / dform[hd]
    return deci

def sep(lat1, lon1, lat2, lon2, hd='d'):
    """
    Calculate seperation between two coordinates in decimal degrees. If using
    longitude in hours set parameter hd to "h".

    Parameters
    ----------
    hd : string, ('h', 'd')
        Hour or degree convention
    """
    dform = {'h': 24 / 360., 'd': 1.}
    lat1, lat2 = _m.radians(lat1), _m.radians(lat2)
    lon1, lon2 = _m.radians(lon1 / dform[hd]), _m.radians(lon2 / dform[hd])
    dlon, dlat = lon1 - lon2, lat1 - lat2
    a = _m.sin(dlat / 2.0)**2 + _m.cos(lat1) * _m.cos(lat2) * _m.sin(dlon / 2.0)**2
    d = 2 * _m.asin(min(1, _m.sqrt(a)))
    return _np.rad2deg(d)

def sep_coords(needle, haystack, min_sep):
    """
    Match a "needle" single (lon, lat) coordinate to a "haystack" list of
    coordinates in decimal degrees. Use sorted lists for best performance.

    Parameters
    ----------
    needle : array like
        List or tuple of (lon, lat) in decimal degrees
    haystack : numpy array
        2 x N list of coordinates in decimal degrees

    Returns
    -------
    sep : numpy array
        Array of seperations compared to original list in radians
    """
    dlon = _np.radians(haystack[:,0]) - _np.radians(needle[0])
    dlat = _np.radians(haystack[:,1]) - _np.radians(needle[1])
    a = _np.square(_np.sin(dlat / 2.0)) + _np.cos(_np.radians(needle[1])) * \
        _np.cos(_np.radians(haystack[:,1])) * _np.square(_np.sin(dlon / 2.0))
    sep = _np.arcsin(_np.minimum(_np.sqrt(a), _np.ones(len(a))))
    return sep

def nearest_match_coords(needle, haystack, min_sep):
    """
    Calculate the nearest source between a "needle" single (lon, lat) coordinate
    and a "haystack" list of coordinates in decimal degrees. Use sorted lists
    for best performance.

    Parameters
    ----------
    needle : array like
        List or tuple of (lon, lat) in decimal degrees
    haystack : numpy array
        2 x N list of coordinates in decimal degrees
    min_sep : number
        Minimum seperation in arcseconds.

    Returns
    -------
    min_index : number
        Array index of nearest object
    min_dist : number
        Distance to nearest matched object
    matchn : number
        Number of matches within the minimum seperation
    """
    min_sep = 2. * _np.pi * min_sep / (360. * 60. * 60.)
    dlon = _np.radians(haystack[:,0]) - _np.radians(needle[0])
    dlat = _np.radians(haystack[:,1]) - _np.radians(needle[1])
    a = _np.square(_np.sin(dlat / 2.0)) + _np.cos(_np.radians(needle[1])) * \
        _np.cos(_np.radians(haystack[:,1])) * _np.square(_np.sin(dlon / 2.0))
    d = _np.arcsin(_np.minimum(_np.sqrt(a), _np.ones(len(a))))
    matchn = len(d[d <= min_sep])
    min_index = d.argmin()
    min_dist = _np.degrees(d.min())
    return matchn, min_index, min_dist

