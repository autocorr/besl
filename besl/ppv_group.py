"""
===========================
PPV Grouping and Clustering
===========================

Functions and routines to perform clustering analysis on the BGPS HCO+/N2H+
molecular line survey.

"""

from __future__ import division
import numpy as np
import pandas as pd
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
    ver = 'v210'
    good_kdars = ['N', 'F', 'T', 'O']

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
        dpdf = read_cat('bgps_kdars_' + self.ver)
        cnums = dpdf[self.ver + 'cnum']
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
        if not self.__scanned:
            raise Exception('Tree has not been built, run `dbscan`.')
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
            if cnum not in df[self.ver + 'cnum'].values:
                self.kdar_skipped += 1
                continue
            ii = df[df[self.ver + 'cnum'] == cnum].index[0]
            cid = self.tree[ii][1]
            cluster_nodes[cid][1].append(kdar)
        # Check unique KDARs
        for cid in cluster_ids:
            if cid == -1:
                continue
            kdar_assoc = cluster_nodes[cid][1]
            kdar_unique = np.unique(kdar_assoc)
            kdar_disagree = ('N' in kdar_unique) & ('F' in kdar_unique)
            outer = 'O' in kdar_unique
            group_f = 0
            if outer:
                group_f = 3
            elif len(kdar_unique) == 1:
                group_f = 1
            elif (len(kdar_unique) > 1) & (not kdar_disagree):
                group_f = 1
            elif (len(kdar_unique) > 1) & kdar_disagree:
                group_f = 2
            cluster_nodes[cid][2] = group_f
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
            cluster_nodes.iteritems() if (v[2] in [1,2,3]) & (k != -1)])
        self.new_kdar_assoc = sum([len(v[0]) - len(v[1]) for k, v in
            cluster_nodes.iteritems() if (v[2] in [1, 2]) & (k != -1)])
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

    def to_df(self):
        """
        Return the cluster-ID and group flag for a BGPS catalog number in
        DataFrame format.
        """
        if not self.__scanned:
            raise Exception('Tree has not been built, run `dbscan`.')
        cols = [self.ver + 'cnum', 'cid', 'group_size', 'group_f',
                'group_kdars']
        table_data = []
        for ix in self.tree.iterkeys():
            cnum = self.df.ix[ix, self.ver + 'cnum']
            cid = self.tree[ix][1]
            group_f = self.cluster_nodes[cid][2]
            if cid == -1:
                group_size = 0
                group_kdars = 0
            else:
                group_size = len(self.cluster_nodes[cid][0])
                group_kdars = len(self.cluster_nodes[cid][1])
            table_data.append([cnum, cid, group_size, group_f, group_kdars])
        self.cluster_df = pd.DataFrame(table_data, columns=cols).sort(cols[0])
        return self.cluster_df


class ClusterRegion(object):
    """
    Create DS9 region files from a `ClusterDBSCAN` instance.

    Parameters
    ----------
    obj : besl.ppv_group.ClusterDBSCAN
    """
    def __init__(self, obj, **kwargs):
        self.obj = obj
        self.angle_rad, self.velo_rad = obj.lims
        self.preamble = 'global color=green font="helvetica 10 normal" ' + \
                        'select=1 highlite=1 edit=1 move=1 delete=1 ' + \
                        'include=1 fixed=0\ngalactic\n'
        self.color_ring = deque(['green', 'red', 'blue', 'yellow', 'magenta',
                                 'cyan', 'orange'])
        self.point_entry = 'x point {l} {b} # text={lcb}{v}:{cid}:{ix}{rcb} color={c}\n'
        self.circle_entry = 'circle {l} {b} ' + str(self.angle_rad) + ' # color={c}\n'
        self.braces = {'lcb': '{', 'rcb': '}'}

    def _write_region(self, out_filen, text):
        with open(out_filen + '.reg', 'w') as f:
            f.write(text)

    def write_all(self):
        """
        Do all write methods with default parameters.
        """
        self.clusters()
        self.clusters(out_filen='ppv_group_nocirc', search_circles=False)
        self.agree_conflict()
        self.single_kdar()
        self.no_kdar()

    def clusters(self, out_filen='ppv_group', search_circles=True):
        """
        Write a DS9 regions file with cluster nodes colored by cluster and circles
        showing the angular search radius, and the nodes velocity visualized.

        Parameters
        ----------
        out_filen : str
            Name of regions files, ends in '.reg' extension.
        search_circles : bool, default False
            Include the circular apertures to extent of angular search radius.

        Attributes
        ----------
        cluster_text : str
            String written to file
        """
        obj = self.obj
        color_ring = self.color_ring
        all_lines = self.preamble
        for cid, params in obj.cluster_nodes.iteritems():
            nodes = params[0]
            kdars = params[1]
            conflict_flag = params[2]
            if cid == -1:
                c = 'black'
            else:
                c = color_ring[0]
            for ii in nodes:
                l = obj.df.ix[ii, 'glon_peak']
                b = obj.df.ix[ii, 'glat_peak']
                v = obj.df.ix[ii, 'all_vlsr']
                cnum = obj.df.ix[ii, self.ver + 'cnum']
                all_lines += self.point_entry.format(l=l, b=b, v=v, cid=cid,
                                                     ix=cnum, c=c, **self.braces)
                if search_circles:
                    all_lines += self.circle_entry.format(l=l, b=b, c=c)
            color_ring.rotate()
        self._write_region(out_filen, text=all_lines)
        self.cluster_text = all_lines

    def agree_conflict(self, out_filen='agree_conflict_group'):
        """
        Write a DS9 regions file with cluster nodes colored by cluster and circles
        showing the angular search radius, and the nodes velocity visualized.

        Parameters
        ----------
        out_filen : str
            Name of regions files, ends in '.reg' extension.

        Attributes
        ----------
        agree_conflict_text : str
            String written to file
        """
        obj = self.obj
        all_lines = self.preamble
        all_colors = {1: 'green', 2: 'red'}
        for cid, params in obj.cluster_nodes.iteritems():
            nodes = params[0]
            kdars = params[1]
            conflict_flag = params[2]
            if (len(kdars) < 2) | (cid == -1):
                continue
            for ii in nodes:
                l = obj.df.ix[ii, 'glon_peak']
                b = obj.df.ix[ii, 'glat_peak']
                c = all_colors[conflict_flag]
                all_lines += self.circle_entry.format(l=l, b=b, c=c)
        self._write_region(out_filen, text=all_lines)
        self.agree_conflict_text = all_lines

    def single_kdar(self, out_filen='single_kdar_group'):
        """
        Write a DS9 regions file for clusters that contain only a single KDAR.

        Parameters
        ----------
        out_filen : str
            Name of regions files, ends in '.reg' extension.

        Attributes
        ----------
        single_kdar_text : str
            String written to file
        """
        obj = self.obj
        all_lines = self.preamble
        for cid, params in obj.cluster_nodes.iteritems():
            nodes = params[0]
            kdars = params[1]
            if len(kdars) == 1:
                for ii in nodes:
                    l = obj.df.ix[ii, 'glon_peak']
                    b = obj.df.ix[ii, 'glat_peak']
                    all_lines += self.circle_entry.format(l=l, b=b, c='grey')
        self._write_region(out_filen, text=all_lines)
        self.single_kdar_text = all_lines

    def no_kdar(self, out_filen='no_kdar_group'):
        """
        Write a DS9 regions file for clusters that contain no KDAR.

        Parameters
        ----------
        out_filen : str
            Name of regions file, ends in '.reg' extension.

        Attributes
        ----------
        no_kdar_text : str
            String written to file
        """
        obj = self.obj
        all_lines = self.preamble
        for cid, params in obj.cluster_nodes.iteritems():
            nodes = params[0]
            kdars = params[1]
            if len(kdars) == 0:
                for ii in nodes:
                    l = obj.df.ix[ii, 'glon_peak']
                    b = obj.df.ix[ii, 'glat_peak']
                    all_lines += self.circle_entry.format(l=l, b=b, c='black')
        self._write_region(out_filen, text=all_lines)
        self.no_kdar_text = all_lines


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

def kdar_flag_region(out_filen='kdar_flags'):
    """
    Write a DS9 region file with the KDA resolution flags at coordinates.

    Parameters
    ----------
    out_filen : str
        Output regions filename with extension '.reg'

    Returns
    -------
    all_lines : str
        Text string for region file content for parsing.
    """
    # Read in data
    df = read_bgps_vel()
    df = df.set_index('v210cnum')
    dpdf = read_cat('bgps_kdars_v210')
    dpdf = dpdf[dpdf['kdar'] != 'U']
    # Preamble
    all_lines = 'global color=green font="helvetica 10 normal" select=1 ' + \
                'highlite=1 edit=1 move=1 delete=1 include=1 fixed=0\n'
    all_lines += 'galactic\n'
    all_colors = {'U': 'black', 'F': 'red', 'N': 'green', 'T': 'orange', 'O':
                  'magenta'}
    # Plot strings
    text_entry = 'text {l} {b} # text={lcb}{flag}{rcb} color={c}\n'
    braces = {'lcb': '{', 'rcb': '}'}
    for cnum, kdar in dpdf.values:
        l = df.ix[cnum, 'glon_peak']
        b = df.ix[cnum, 'glat_peak'] - 0.005
        c = all_colors[kdar]
        all_lines += text_entry.format(l=l, b=b, c=c, flag=kdar, **braces)
    with open(out_filen + '.reg', 'w') as f:
        f.write(all_lines)
    return all_lines


