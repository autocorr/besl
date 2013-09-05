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
from .catalog import read_bgps_vel, read_dpdf


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

class ClusterDBSCAN(object):
    """
    Implimentation of DBSCAN cluster recognition algorithm.

    Parameters
    ----------
    df : pd.DataFrame
    cols : list, default ['glon_peak', 'glat_peak', 'vlsr']
        List of column names for Galactic longitude, latitude, and velocity in
        decimal degrees and km/s.
    flag_col : str, 'vlsr_'
        Flag column to select good velocities from.
    lims : list, default [0.1, 3.5]
        Coordinates search radius for [angle, velocity] where angle is written
        in decimal degrees.
    min_points : number, default 1
        Minimum number of points in a nodes neighborhood in order to associate
        it with a cluster.

    Returns
    -------
    tree : dict
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
        # Initialize tree and BGPS
        self.good_kdars = ['N', 'F', 'T', 'O']
        self.cluster_id = 0
        self.read_data()
        self.velo_to_ang = lims[0] / lims[1]
        self.__scanned = False

    def dbscan(self):
        df = self.df
        for ii in df.index:
            # Mark as visited
            self.tree[ii][0] = True
            neighbors = self.region_query(ii)
            self.tree[ii][2].extend(neighbors)
            if len(neighbors) <= self.min_points:
                self.tree[ii][1] = -1
            else:
                self.cluster_id += 1
                self.expand_cluster(ii, neighbors)
        self.__scanned = True

    def expand_cluster(self, ix, neighbors):
        self.tree[ix][1] = self.cluster_id
        neighbors = set(neighbors)
        visited_neighbors = set()
        # Recursively search neighbors
        import ipdb; ipdb.set_trace()
        while neighbors:
            ii = neighbors.pop()
            visited = self.tree[ii][0]
            node_id = self.tree[ii][1]
            visited_neighbors.add(ii)
            if not visited:
                # Mark as visited
                self.tree[ii][0] = True
                branch = set(self.region_query(ii))
                # Union branch to current set of neighbors if not visited
                if len(branch) > self.min_points:
                    neighbors.update(branch.difference(visited_neighbors))
            # If not yet a member, assign to cluster
            if node_id == 0:
                self.tree[ii][1] = self.cluster_id

    def region_query(self, ix):
        df = self.df
        lim_a, lim_v = self.lims
        cols = self.cols
        l, b, v = df.ix[ix, cols].values
        # Select box
        df = df[(df[cols[0]] > l - lim_a) &
                (df[cols[0]] < l + lim_a) &
                (df[cols[1]] > b - lim_a) &
                (df[cols[1]] < b + lim_a) &
                (df[cols[2]] > v - lim_v) &
                (df[cols[2]] < v + lim_v)]
        # Coordinate query
        # FIXME doesn't do correct latitude transformation
        neighbors = df[(df[cols[0]] - l)**2 + (df[cols[1]] - b)**2 +
                       self.velo_to_ang**2 * (df[cols[2]] - v)**2 <
                       lim_a**2].index.values
        return neighbors

    def read_data(self):
        df = read_bgps_vel()
        df = df[(df[self.cols[2]].notnull()) & (df[self.flag_col] > 0)]
        dpdf = read_dpdf()
        kdars = dpdf[6].data['KDAR']
        cnums = dpdf[1].data['CNUM']
        # Check well resolved KDAs
        self.good_cnums = cnums[np.in1d(kdars, self.good_kdars)]
        self.good_cnums_kdars = kdars[np.in1d(kdars, self.good_kdars)]
        # Assign as instance variables
        self.dpdf = dpdf
        self.df = df
        self.tree = {ix : [False, 0, []] for ix in df.index}

    def analysis(self):
        if self.__scanned:
            tree = self.tree
            df = self.df
            # Number of clusters
            cluster_ids = np.unique([row[1] for row in tree.values()])
            n_clusters = cluster_ids.shape[0]
            # Cluster nodes
            cluster_nodes = {}
            for cid in cluster_ids:
                nodes = []
                for ix in df.ix:
                    if tree[ix][1] == cid:
                        nodes.extend(tree[ix][2])
                cluster_nodes[cid] = list(np.unique(nodes))
            # TODO
            # KDAR nodes in cluster
            all_core_nodes = np.ravel(cluster_nodes.values())
            good_cnums = self.good_cnums
            good_cnums_kdars = self.good_cnums_kdars
            #kdar_cnums = 
            # Number of nodes in cluster with KDARs
            # Number of nodes in clusters with KDAR conflicts
            # Assign and save results
            results = {}
            results['ids'] = cluster_ids
            results['n_clusters'] = n_clusters
            results['nodes'] = cluster_nodes
            self.results = results
        else:
            raise Exception('Tree has not been built, run `dbscan`.')


