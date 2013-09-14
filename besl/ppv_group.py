"""
===========================
PPV Grouping and Clustering
===========================

Functions and routines to perform clustering analysis on the BGPS HCO+/N2H+
molecular line survey.

"""

from __future__ import division
import numpy as np
import cPickle as pickle
from collections import deque
from multiprocessing import Pool
from rtree import index
from .catalog import read_bgps_vel, read_cat


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
    cols : list, default ['glon_peak', 'glat_peak', 'all_vlsr']
        List of column names for Galactic longitude, latitude, and velocity in
        decimal degrees and km/s.
    flag_col : str, 'vlsr_f'
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

    Attributes
    ----------
    kdar_conflict_nodes : number
    kdar_conflict_clusters : number
    kdar_agree_nodes : number
    kdar_agree_clusters : number
    kdar_span_nodes : number
    new_kdar_assoc : number
    conflict_frac : number
    cluster_ids : list
    n_clusters : number
    n_core_nodes : number
    self.cluster_nodes : dict
    """
    def __init__(self,
                 cols=['glon_peak', 'glat_peak', 'all_vlsr'],
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
        """
        Group clusters of points in PPV-space using the `DBSCAN` algorithm.
        """
        df = self.df
        for ii in df.index:
            if self.tree[ii][0]:
                continue
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
        """
        Recursively search the nodes of `tree` if they are flagged as
        univisited and add them to the current `cluster_id`. The `self.tree` is
        modified in place.

        Parameters
        ----------
        ix : number
            Node dataframe index number
        neighbors : array-like
            Node neighbor indices
        """
        self.tree[ix][1] = self.cluster_id
        neighbors = set(neighbors)
        visited_neighbors = set()
        # Recursively search neighbors
        while neighbors:
            ii = neighbors.pop()
            visited = self.tree[ii][0]
            node_id = self.tree[ii][1]
            visited_neighbors.add(ii)
            if not visited:
                # Mark as visited
                self.tree[ii][0] = True
                branch = set(self.region_query(ii))
                self.tree[ii][2].extend(branch)
                # Union branch to current set of neighbors if not visited
                if len(branch) > self.min_points:
                    neighbors.update(branch.difference(visited_neighbors))
            # If not yet a member, assign to cluster
            if node_id == 0:
                self.tree[ii][1] = self.cluster_id

    def region_query(self, ix):
        """
        Search the node's coordinate neighborhood of index `ix` and return
        indices of neighbors.

        Parameters
        ----------
        ix : number
            Node dataframe index number

        Returns
        -------
        neighbors : np.array
            Array of neighbor node indices
        """
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
        """
        Read and process data for DPDFs, velocity catalog, and select sources
        with well-resolved KDAs.
        """
        df = read_bgps_vel()
        df = df[df[self.cols[2]].notnull()]
        dpdf = read_cat('bgps_kdars_v210')
        cnums = dpdf['v210cnum']
        kdars = dpdf['kdar']
        # Check well resolved KDAs
        dpdf_mask = np.in1d(kdars, self.good_kdars)
        # Assign as instance variables
        self.good_cnums = zip(cnums[dpdf_mask], kdars[dpdf_mask])
        self.dpdf = dpdf
        self.df = df
        self.tree = {ix : [False, 0, []] for ix in df.index}

    def analysis(self, verbose=False):
        """
        Analyze the built tree. Requires `dbscan` to have been already run.

        Parameters
        ----------
        verbose : bool
            Print results to terminal.
        """
        if self.__scanned:
            tree = self.tree
            df = self.df
            # Number of clusters
            cluster_ids = np.unique([row[1] for row in tree.values()])
            n_clusters = cluster_ids.shape[0]
            # Cluster nodes
            cluster_nodes = {ix : [[], [], 0] for ix in cluster_ids}
            for cid in cluster_ids:
                nodes = []
                for ix in df.index:
                    if tree[ix][1] == cid:
                        nodes.extend(tree[ix][2])
                cluster_nodes[cid][0].extend(np.unique(nodes))
            # Nodes in clusters
            core_nodes = cluster_nodes.copy()
            del core_nodes[-1]
            n_core_nodes = np.ravel(core_nodes.values()).shape[0]
            # KDAR nodes in cluster
            good_cnums = self.good_cnums
            self.kdar_skipped = 0
            for cnum, kdar in good_cnums:
                # Should only continue if more stringent velocity flags
                if cnum not in df['v210cnum'].values:
                    self.kdar_skipped += 1
                    continue
                ii = df[df['v210cnum'] == cnum].index[0]
                cid = self.tree[ii][1]
                cluster_nodes[cid][1].append(kdar)
            # Check unique KDARs
            for cid in cluster_ids:
                kdar_assoc = cluster_nodes[cid][1]
                kdar_confused = np.unique(kdar_assoc).shape[0]
                if kdar_confused == 1:
                    cluster_nodes[cid][2] = 1
                elif kdar_confused > 1:
                    cluster_nodes[cid][2] = 2
            # Number of nodes in clusters with KDAR conflicts
            self.kdar_conflict_nodes = sum([len(v[0]) for k, v in
                cluster_nodes.iteritems() if (v[2] == 2) & (k != -1)])
            self.kdar_conflict_clusters = len([len(v[0]) for k, v in
                cluster_nodes.iteritems() if (v[2] == 2) & (k != -1)])
            # Number of nodes in clusters with agreeing KDARs
            self.kdar_agree_nodes = sum([len(v[0]) for k, v in
                cluster_nodes.iteritems() if (v[2] == 1) & (len(v[1]))])
            self.kdar_agree_clusters = len([len(v[0]) for k, v in
                cluster_nodes.iteritems() if (v[2] == 1) & (len(v[1]))])
            # Number of nodes in cluster with KDARs
            self.kdar_span_nodes = sum([len(v[0]) for k, v in
                cluster_nodes.iteritems() if (v[2] in [1,2]) & (k != -1)])
            self.new_kdar_assoc = self.kdar_span_nodes + \
                len(cluster_nodes[-1][1]) - np.shape(good_cnums)[0] + \
                self.kdar_skipped
            self.conflict_frac = self.kdar_conflict_nodes / \
                (self.kdar_agree_nodes + self.kdar_conflict_nodes)
            # Assign and save results
            self.cluster_ids = cluster_ids
            self.n_clusters = n_clusters
            self.n_core_nodes = n_core_nodes
            self.cluster_nodes = cluster_nodes
            if verbose:
                print """
                -- Results:
                {0:<10} : Clusters
                {1:<10} : Cluster (core) nodes
                {2:<10} : Nodes in clusters containing KDAR clumps
                {3:<10} : Net new KDAR associations for nodes in clusters
                {4:<10} : Nodes in clusters containing KDAR conflicts
                {5:<10} : Nodes in clusters containing KDAR agreements
                {6:<10f} : Ratio of conflict nodes to all multi-KDAR nodes
                """.format(n_clusters, n_core_nodes, self.kdar_span_nodes,
                           self.new_kdar_assoc, self.kdar_conflict_nodes,
                           self.kdar_agree_nodes, self.conflict_frac)
        else:
            raise Exception('Tree has not been built, run `dbscan`.')

def grid_calc(lims=[0.05, 0.2, 1, 4], points=[10, 10], out_filen='obj_grid'):
    """
    Run DBSCAN for a grid of angle and velocity search distances. Uses
    multiprocssesing by default.

    Parameters
    ----------
    lims : list
        List of limits for search [angle_min, angle_max, velocity_min,
        velocity_max].
    points : list
        Number of grid points to sample, end-inclusive.
    out_filen : str
        Filename for pickled object grid, ends in '.pickle' extension.

    Returns
    -------
    obj_grid : np.array
        Object array of `ClusterDBSCAN` instances
    X, Y : np.array
        Result of np.meshgrid over arrays of sample points
    """
    assert (len(lims) == 4) & (len(points) == 2)
    assert (points[0] >= 2) & (points[1] >= 2)
    x = np.linspace(lims[0], lims[1], points[0])
    y = np.linspace(lims[2], lims[3], points[1])
    X, Y = np.meshgrid(x, y)
    limits = np.dstack([X, Y]).reshape(-1, 2)
    clusters = (ClusterDBSCAN(lims=l) for l in limits)
    # Compute clusters with multiprocessing
    pool = Pool(processes=6)
    obj_grid = pool.imap(wrapper, clusters)
    pool.close()
    pool.join()
    # Reshape to grid
    obj_grid = np.reshape([obj for obj in obj_grid], X.shape)
    with open(out_filen + '.pickle', 'wb') as f:
        pickle.dump([obj_grid, X, Y], f)
    return obj_grid, X, Y

def wrapper(c):
    """
    Wrapper function on the top-level domain so that the object is pickleable
    for `multiprocessing.Pool`.
    """
    c.dbscan()
    c.analysis(verbose=True)
    return c

def reduce_obj_grid(filen='obj_grid', out_filen='obj_props'):
    """
    Extract parameters from object grid into a dictionary of matricies for easy
    plotting and analysis.

    Parameters
    ----------
    filen : str, default 'obj_grid'
        Filename of the `obj_grid` pickle ending in the '.pickle' extension.
    out_filen : str, default 'obj_props'
        Filename of the reduced object parameter dictionary, ends in the
        '.pickle' extension.

    Returns
    -------
    obj_dict : dict
        Dictionary of the analysis values over each tree in `obj_grid`.
    """
    obj_grid, X, Y = pickle.load(open(filen + '.pickle', 'rb'))
    obj_dict = {}
    obj_dict['angle'] = X
    obj_dict['velo'] = Y
    props = [('n_clusters', lambda c: c.n_clusters),
             ('kdar_conflict_nodes', lambda c: c.kdar_conflict_nodes),
             ('kdar_conflict_clusters', lambda c: c.kdar_conflict_clusters),
             ('kdar_agree_nodes', lambda c: c.kdar_agree_nodes),
             ('kdar_agree_clusters', lambda c: c.kdar_agree_clusters),
             ('kdar_span_nodes', lambda c: c.kdar_span_nodes),
             ('new_kdar_assoc', lambda c: c.new_kdar_assoc),
             ('conflict_frac', lambda c: c.conflict_frac)]
    for key, method in props:
        obj_dict[key] = np.reshape(map(method, obj_grid.flat), X.shape)
    with open(out_filen + '.pickle', 'wb') as f:
        pickle.dump(obj_dict, f)

def cluster_region(obj, out_filen='ppv_group'):
    """
    Write a DS9 regions file with cluster nodes colored by cluster and circles
    showing the angular search radius, and the nodes velocity visualized.

    Parameters
    ----------
    obj : besl.ppv_group.ClusterDBSCAN
        ClusterDBSCAN instance to extract parameters from
    out_filen : str
        Name of regions files, ends in '.reg' extension.

    Returns
    -------
    all_lines : str
        String written to file
    """
    angle_rad, velo_rad = obj.lims
    all_lines = 'global color=green font="helvetica 10 normal" select=1 ' + \
                'highlite=1 edit=1 move=1 delete=1 include=1 fixed=0\n'
    all_lines += 'galactic\n'
    colors = deque(['green', 'red', 'yellow', 'blue', 'magenta', 'cyan'])
    point_entry = 'x point {l} {b} # text={lcb}{v},{cid}{rcb} color={c}\n'
    circle_entry = 'circle {l} {b} ' + str(angle_rad) + ' # color={c}\n'
    braces = {'lcb': '{', 'rcb': '}'}
    for cid, params in obj.cluster_nodes.iteritems():
        c = colors[0]
        nodes = params[0]
        kdars = params[1]
        conflict_flag = params[2]
        for ii in nodes:
            l = obj.df.ix[ii, 'glon_peak']
            b = obj.df.ix[ii, 'glat_peak']
            v = obj.df.ix[ii, 'all_vlsr']
            all_lines += point_entry.format(l=l, b=b, v=v, cid=ii, c=c,
                                            **braces)
            all_lines += circle_entry.format(l=l, b=b, c=c)
        colors.rotate()
    with open(out_filen + '.reg', 'w') as f:
        f.write(all_lines)
    return all_lines


