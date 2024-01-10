from typing import Union
from itertools import combinations, product

import numpy as np


class Peak():
    def __init__(self, middle: int, seq_len: int, window_size: int) -> "Peak":
        '''
        Peak object for keeping track of points with a window (i.e. border/buffer) surrounding them.
        Expected behaviour is that `self.middle` is always in the middle position between `self.five_side` and `self.three_side`.
        - `five_side`  : 5' side of the window.
        - `three_side` : 3' side of the window.

        Alternate constructors include: `from_calc_middle` and `from_edges`.
        '''
        self.middle = middle
        self.seq_len = seq_len
        self.window_size = window_size
        self.split = False

        self.z_score: int = 0
        self.g_score: int = 0
        self.d_score: int = 0
        self.decision: bool = None

        five_side = self.middle - self.window_size // 2
        three_side = self.middle + self.window_size // 2

        if five_side < 0:
            self.split = True
            five_side += self.seq_len-1
        if three_side > self.seq_len-1:
            self.split = True
            three_side -= self.seq_len-1

        self.five_side = five_side
        self.three_side = three_side

    @classmethod
    def from_calc_middle(cls, five_side, three_side, seq_len, window_size) -> "Peak":
        '''Calculates the middle of five_side and three_side and uses that as the middle'''
        return cls(cls.calc_middle(five_side, three_side, seq_len), seq_len, window_size)

    @classmethod
    def from_edges(cls, five_side, three_side, seq_len) -> "Peak":
        '''
        Calculates the middle of five_side and three_side and uses that as the middle 
        AND uses these edges to determine the window_size.
        '''
        return cls(cls.calc_middle(five_side, three_side, seq_len), seq_len, cls.calc_dist_points(five_side, three_side, seq_len))

    @staticmethod
    def calc_dist_points(a: int, b: int, curve_size: int) -> int:
        '''Calculate the distance of a to b in bp on a circular chromosome.'''
        dist_1 = max(a, b) - min(a, b)
        dist_2 = min(a, b) + curve_size - max(a, b)
        return min(dist_1, dist_2)

    @staticmethod
    def get_peaks_to_merge(peaks: list["Peak"]) -> list[tuple["Peak", "Peak"]]:
        """
        Get the indices that give the same value in the curve and are in eachothers window.
        Input:
        - `peaks`          : list of `Peaks` of a curve
        Return:
        - `peaks_to_merge` : Nested list. Each sublist contains the two Peaks that have to be merged
        """
        peaks_to_merge = []
        for peak_i, peak_j in combinations(peaks, 2):
            if peak_i.contains_point(peak_j.middle):
                peaks_to_merge.append((peak_i, peak_j))
        return peaks_to_merge

    @staticmethod
    def match_peaks(peaks_x: list["Peak"], peaks_y: list["Peak"]) -> list[tuple["Peak", "Peak"]]:
        '''Nested list of peaks from x that line up with peaks from y.'''
        matched_peaks = []
        for peak_x, peak_y in product(peaks_x, peaks_y):
            if peak_x.intersecting_windows(peak_y):
                matched_peaks.append((peak_x, peak_y))
        return matched_peaks

    @staticmethod
    def get_adjacency_matrix(peaks_a: list[Union["Peak", int]], peaks_b: list[Union["Peak", int]] = None, seq_len: int = None) -> np.ndarray:
        """
        Gets adjacency matrix for given list of `Peak` or `int` objects.
        The matrix can be between a list and the elements of itself or between the elements of two lists.
        All elements in each list must be of the same type.
        Elements in `peaks_a` must be the same type as those in `peaks_b`.
        If elements are integers, `seq_len` must be provided.
        """

        # Error handling and variable initialisation
        are_integers = True if isinstance(peaks_a[0], int) else False
        if are_integers and seq_len is None:
            raise ValueError('Provided list of integers, but did not provide `seq_len`.')
        if not all(isinstance(x, type(peaks_a[0])) for x in peaks_a[1:]):
            raise ValueError('Not all elements in `peaks_a` are of the same type.')

        if peaks_b is not None:
            if not all(isinstance(x, type(peaks_b[0])) for x in peaks_b[1:]):
                raise ValueError('Not all elements in `peaks_b` are of the same type.')
            if not isinstance(peaks_b[0], type(peaks_a[0])):
                raise ValueError('Elements is `peaks_a` are not the same type as `peaks_b`.')

            adj_mat = np.zeros((len(peaks_a), len(peaks_b)))
            iterator = product(enumerate(peaks_a), enumerate(peaks_b))
        else:
            adj_mat = np.zeros((len(peaks_a), len(peaks_a)))
            iterator = combinations(enumerate(peaks_a), r=2)

        for (i_a, a), (i_b, b) in iterator:
            dist = Peak.calc_dist_points(a, b, seq_len) if are_integers else Peak.calc_dist_points(
                a.middle, b.middle, a.seq_len)
            adj_mat[i_a, i_b] = dist
            if peaks_b is None:
                adj_mat[i_b, i_a] = dist
        return adj_mat

    @staticmethod
    def get_connected_components_idx(vertices: list[int], adj_mat: np.ndarray, threshold: int) -> list:
        """
        Recursively find connected groups in an undirected graph. 

        Parameters:
        - `vertices`  : list of indexes of whatever you want to group
        - `adj_mat`   : adjacency matrix of len(vertices) x len(vertices).
        - `threshold` : maximum size of an edge before two vertices are considered connected.
        """

        def _DFS_recurse(idx: int, adj_mat: np.ndarray, visited: list[bool], connected_list: list[int], threshold: int):
            """Private: used by get_connected_groups_idx for recursion"""
            visited[idx] = True
            connected_list.append(idx)
            for i in range(len(visited)):
                if adj_mat[i][idx] < threshold and not visited[i]:
                    _, _, visited, connected_list, _ = _DFS_recurse(i, adj_mat, visited, connected_list, threshold)
            return idx, adj_mat, visited, connected_list, threshold

        visited = [False] * len(vertices)
        connected_groups_idx = []
        for i in range(len(vertices)):
            if not visited[i]:
                group = []
                _, _, visited, group, _ = _DFS_recurse(i, adj_mat, visited, group, threshold=threshold)
                connected_groups_idx.append(group)
        return connected_groups_idx

    @staticmethod
    def select_connected_groups(peaks: list["Peak"], adj_mat: np.ndarray, threshold: int) -> list[list["Peak"]]:
        """
        Gets the connected components in an undirected graph and limits the within-distance of clusters to < threshold*3.

        Parameters:
        - `peaks`     : list of Peak objects you want to group
        - `adj_mat`   : adjacency matrix for each Peak to each other Peak
        - `threshold` : maximum distance between peaks on the sequence for which they are considered connected.
        """
        cg_i = Peak.get_connected_components_idx([i for i in range(len(peaks))], adj_mat, threshold)
        groups_to_process = [(group, threshold) for group in cg_i]
        accepted_groups_idx = []
        flag = False
        while len(groups_to_process) != 0:
            curr_group, curr_threshold = groups_to_process.pop(0)
            # check within distance of group
            for i, j in combinations(curr_group, r=2):
                # NOTE: Tune the threshold!!!!! especially check for large and realistic sequences!
                if adj_mat[i][j] > threshold*1.5:
                    flag = True
                    break
            # got flagged as too wide.
            if flag and curr_threshold > 0.25*threshold:
                # get peaks that correspond to the indexes in the flagged group
                curr_peaks = [peaks[i] for i in curr_group]
                mat = Peak.get_adjacency_matrix(curr_peaks)
                new_groups = Peak.get_connected_components_idx(curr_peaks, mat, curr_threshold)
                # revert indexing to that of the peaks list
                new_groups = [[curr_group[i] for i in group] for group in new_groups]
                # add new groups for processing
                groups_to_process.extend([(new_group, 0.75*curr_threshold) for new_group in new_groups])
                flag = False
            else:
                accepted_groups_idx.append(curr_group)

        connected_groups_vals = [[peaks[i] for i in idx_group] for idx_group in accepted_groups_idx]
        return connected_groups_vals

    @staticmethod
    def calc_middle(a: Union[int, "Peak"], b: Union[int, "Peak"], curve_size: int = None) -> int:
        """
        Calculate the distance between Peak a and b on circular DNA.
        Returns a new Peak.middle (int) in the middle of a and b
        """

        seq_len = a.seq_len if curve_size is None and isinstance(a, Peak) else curve_size
        seq_len = b.seq_len if seq_len is None and isinstance(b, Peak) else seq_len

        a, b = int(a), int(b)

        if seq_len is None:
            raise ValueError("No curve size found in any parameters. This method is for circular DNA only.")

        dist_1 = max(a, b) - min(a, b)
        dist_2 = min(a, b) + seq_len-1 - max(a, b)

        if dist_2 < dist_1:
            merged_idx = min(a, b) - dist_2//2
            if merged_idx < 0:
                merged_idx = seq_len-1 + merged_idx
        else:
            merged_idx = min(a, b) + dist_1//2
        return merged_idx

    def intersecting_windows(self, other: "Peak") -> bool:
        '''T|F wether self's window intersects with other's window.'''
        if not isinstance(other, Peak):
            raise ValueError(f'other is a {type(other)}, must be a Peak object')
        x = Peak.calc_dist_points(self.middle, other.middle, self.seq_len)
        y = self.window_size//2 + other.window_size//2
        return x < y

    def contains_point(self, point: Union["Peak", int, float]) -> bool:
        '''T|F wether a point is within a Peak's window'''
        p = point.middle if isinstance(point, Peak) else point
        return self.window_size // 2 >= Peak.calc_dist_points(self.middle, p, self.seq_len)

    def calc_dist(self, point_b: Union["Peak", int]):
        """Calculate the distance between self and another point on a cicular chromosome."""
        b = point_b.middle if isinstance(point_b, Peak) else point_b
        return Peak.calc_dist_points(self.middle, b, self.seq_len)

    def __hash__(self) -> int:
        '''very simple hash function. Does not include scores'''
        return hash(
            (self.middle, self.seq_len, self.window_size, self.split,
             self.five_side, self.three_side)
        )

    def __eq__(self, other: "Peak") -> bool:
        '''Simple equality check. Is not used in __lt__ and __gt__'''
        if not isinstance(other, Peak):
            return False
        return self.middle == other.middle \
            and self.seq_len == other.seq_len \
            and self.window_size == other.window_size \
            and self.split == other.split \
            and self.five_side == other.five_side \
            and self.three_side == other.three_side

    def __contains__(self, point: Union["Peak", int, float]) -> bool:
        """Return point in self's window. See: contains_point"""
        return self.contains_point(point)

    def __repr__(self) -> str:
        """Return repr(self)"""
        return f'Peak(middle={self.middle}, window_size={self.window_size}, seq_len={self.seq_len})'

    def __str__(self) -> str:
        """Return str(self)"""
        return str(self.middle)

    def __neg__(self) -> "Peak":
        """Returns -self with a Peak with negation of self.middle"""
        return Peak(-self.middle, self.seq_len, self.window_size)

    def __lt__(self, other: "Peak") -> bool:
        """Return self < other"""
        return self.middle < other.middle

    def __gt__(self, other: "Peak") -> bool:
        """Return self > other"""
        return self.middle > other.middle

    def __add__(self, other: Union["Peak", int, float]) -> Union["Peak", int, float]:
        """Return self + other"""
        if isinstance(other, Peak):
            return Peak(self.middle + other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle + other

    def __sub__(self, other: Union["Peak", int, float]) -> Union["Peak", int, float]:
        """Return self - other"""
        if isinstance(other, Peak):
            return Peak(self.middle - other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle - other

    def __radd__(self, other: Union["Peak", int, float]) -> Union["Peak", int, float]:
        """Return other + self"""
        return self.__add__(other)

    def __rsub__(self, other: Union["Peak", int, float]) -> Union["Peak", int, float]:
        """Return other - self"""
        return -self + other

    def __int__(self):
        """Return int of self.middle"""
        return int(self.middle)

    def __float__(self):
        """Return float of self.middle"""
        return float(self.middle)
