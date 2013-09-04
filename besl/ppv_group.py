"""
===========================
PPV Grouping and Clustering
===========================

Functions and routines to perform clustering analysis on the BGPS HCO+/N2H+
molecular line survey.

"""

from __future__ import division
import numpy as np
from rtree import index
from .catalog import read_bgps_vel


def build_tree(df=None, cols=['glon_peak','glat_peak','vlsr'],
        flag_col='vlsr_f', box_size=[0.2, 0.2, 6]):
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
    box_size : list
        Query region for each data point, specifies the halfwidth of the box.

    Returns
    -------
    ix3d : rtree.index.Index
    """
    # Read in BGPS
    if df is None:
        df = read_bgps_vel()
    if flag_col is None:
        df = df[(df[cols[2]].notnull())]
    else:
        df = df[(df[cols[2]].notnull()) & (df[flag_col] > 0)]
    # Index properties
    props = index.Property()
    props.dimension = 3
    ix3d = index.Index(properties=props)
    # Add coordinates
    box_size = np.array(box_size)
    for ii in df.index:
        coords = df.ix[ii, cols].values
        box = np.ravel([coords - box_size, coords + box_size])
        ix3d.insert(ii, box)
    return ix3d

def dbscan(df=None, cols=['glon_peak','glat_peak','vlsr'],
           flag_col='vlsr_f', lims=[0.1, 3.5]):
    """
    Implimentation of DBSCAN cluster recognition algorithm.

    Parameters
    ----------
    df : pd.DataFrame
    cols : list
        List of column names for Galactic longitude, latitude, and velocity in
        decimal degrees and km/s.
    flag_col : str
        Flag column to select good velocities from.
    lims : list
        Coordinates search radius for [angle, velocity] where angle is written
        in decimal degrees.
    min_points : number
        Minimum number of points in a nodes neighborhood in order to associate
        it with a cluster.

    Returns
    -------
    tree : np.recarray
    """
    # Read in BGPS
    if df is None:
        df = read_bgps_vel()
    if flag_col is None:
        df = df[(df[cols[2]].notnull())]
    else:
        df = df[(df[cols[2]].notnull()) & (df[flag_col] > 0)]
    # Construct tree
    tree = {ix : [0, 0, []] for ix in df.index}
    cluster_id = 0
    for ii in df.index:
        if len(neighbors) < min_points:
            pass
        else:
            cluster_id += 1
            expand_cluster(ii, neighbors, cluster_id, eps, min_points)
    return

class ClusterDBSCAN(object):
    """
    """
    def __init__(self,
                 cols=['glon_peak', 'glat_peak', 'vlsr'],
                 flag_col='vlsr_f',
                 lims=[0.1, 3.5],
                 min_points=1,
                 **kwargs):
        # Set values
        self.cols = cols
        self.flag_col = flag_col
        self.lims = lims
        self.min_points = 1
        self.cluster_id = 0
        # Initialize tree and BGPS
        self.read_data()
        self.eps = lims[0] / lims[1]

    def dbscan(self):
        df = self.df
        for ii in df.index:
            neighbors = self.reqion_query(coord)
            if len(neighbors) < min_points:
                pass
            else:
                self.cluster_id += 1
                expand_cluster(ii, neighbors, min_points)

    def expand_cluster(self):
        pass

    def region_query(self, coord):
        df = self.df
        angle, velo = self.lims
        cols = self.cols
        # Select box
        df = df[(df[cols[0]] > coord[0] - angle) &
                (df[cols[0]] < coord[0] + angle) &
                (df[cols[1]] > coord[1] - angle) &
                (df[cols[1]] < coord[1] + angle) &
                (df[cols[2]] > coord[2] - velo) &
                (df[cols[2]] < coord[2] + velo)]
        # Coordinate query
        return neighbors

    def read_data(self):
        df = read_bgps_vel()
        self.df = df[(df[self.cols[2]]) & (df[self.flag_col] > 0)]
        self.tree = {ix : [0, 0, []] for ix in df.index}






