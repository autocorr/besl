"""
===========================
PPV Grouping and Clustering
===========================

Functions and routines to perform clustering analysis on the BGPS HCO+/N2H+
molecular line survey.

"""

from rtree import index
from .catalog import read_bgps_vel


def build_tree(df=None, cols=['glon_peak','glat_peak','vlsr'],
        flag_col='vlsr_f'):
    """
    Build the R-tree datastructure for the PPV values. If the function is not
    passed any arguments, then the BGPS HCO+/N2H+ catalog is used.

    Parameters
    ----------
    df : pd.DataFrame
    cols : list
        List of column names for Galactic longitude, latitude, and velocity in
        decimal degrees and km/s.
    flag_col : str
        Flag column to select good velocities from.

    Returns
    -------
    ix3d : rtree.
    """
    # Read in BGPS
    if df is None:
        df = read_bgps_vel()
    df = df[(bgpv[cols[2]].notnull()) & (bgpv[flag_col] > 0)]
    # Create index
    props = index.Property()
    props.dimension = 3
    ix3d = index.Index(properties=props)
    # Add coordinates
    for ii in df.index:
        ix3d.add(ii, df.ix[ii, cols].values)
    return ix3d


