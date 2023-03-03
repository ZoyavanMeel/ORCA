from typing import Union
from itertools import combinations, product

import numpy as np

class Peak():
    def __init__(self, middle: int, seq_len: int, window_size: int):
        '''
        five_side  : 5' side of the window
        three_side : 3' side of the window.
        '''
        self.middle = middle
        self.seq_len = seq_len
        self.window_size = window_size
        self.split = False

        self.z_score = 0
        self.g_score = 0
        self.d_score = 0
        self.decision = 0

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
    def from_calc_middle(cls, five_side, three_side, seq_len, window_size):
        '''Calculates the middle of five_side and three_side and uses that as the middle'''
        return cls( cls.get_middle(five_side, three_side, seq_len), seq_len, window_size )
    

    @classmethod
    def from_edges(cls, five_side, three_side, seq_len):
        '''
        Calculates the middle of five_side and three_side and uses that as the middle 
        AND uses these edges to determine the window_size.
        '''
        return cls( cls.get_middle(five_side, three_side, seq_len), seq_len, cls.calc_dist(five_side, three_side, seq_len) )


    @staticmethod
    def calc_dist( a: int, b: int, curve_size: int) -> int:
        '''Calculate the distance of self to other in bp'''
        dist_1 = max(a, b) - min(a, b)
        dist_2 = min(a, b) + curve_size - max(a, b)
        return min(dist_1, dist_2)


    @staticmethod
    def get_peaks_to_merge(peaks: list["Peak"]) -> list[tuple["Peak", "Peak"]]:
        """
        Get the indeces that give the same value in the curve and are in eachothers window.
        Input:
        - `peaks`          : list of `Peaks` of a curve
        Return:
        - `peaks_to_merge` : Nested list. Each sublist contains the two Peaks that have to be merged
        """
        peaks_to_merge = []
        for peak_i, peak_j in combinations(peaks, 2):
            if peak_i.contains_point(peak_j.middle):
                peaks_to_merge.append( (peak_i, peak_j) )
        return peaks_to_merge


    @staticmethod
    def match_peaks(peaks_x: list["Peak"], peaks_y: list["Peak"]) -> list[tuple["Peak", "Peak"]]:
        '''Nested list of peaks from x that line up with peaks from y.'''
        matched_peaks = []
        for peak_x, peak_y in product(peaks_x, peaks_y):
            if peak_x.intersecting_windows(peak_y):
                matched_peaks.append( (peak_x, peak_y) )
        return matched_peaks


    @staticmethod
    def get_adjacency_matrix(peaks_a: list["Peak"], peaks_b: list["Peak"] = None, seq_len: int = None) -> np.ndarray:
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
            if not all(isinstance(x, type(peaks_a[0])) for x in peaks_a[1:]):
                raise ValueError('Not all elements in `peaks_b` are of the same type.')
            if not isinstance(peaks_b[0], type(peaks_a[0])):
                raise ValueError('Elements is `peaks_a` are not the same type as `peaks_b`.')

            adj_mat = np.zeros((len(peaks_a), len(peaks_b)))
            iterator = product(enumerate(peaks_a), enumerate(peaks_b))
        else:
            adj_mat = np.zeros((len(peaks_a), len(peaks_a)))
            iterator = combinations(enumerate(peaks_a), r=2)

        # The function
        for (i_a, a), (i_b, b) in iterator:
            dist = Peak.calc_dist(a, b, seq_len) if are_integers else Peak.calc_dist(a.middle, b.middle, a.seq_len)
            adj_mat[i_a, i_b] = dist
            if peaks_b is None:
                adj_mat[i_b, i_a] = dist
        return adj_mat


    @staticmethod
    def get_connected_groups(peaks: list, adj_mat: np.ndarray, threshold: int) -> list:
        """Recursively find connected groups in an undirected graph"""

        def get_connected_groups_init(peaks: list, adj_mat: np.ndarray, threshold: int) -> list:
            """Private: Groups initial indices of `peaks`."""
            visited = [False] * len(peaks)
            connected_groups_idx = []
            for i in range(len(peaks)):
                if not visited[i]:
                    group = []
                    _, _, visited, group, _ = DFS_recurse(i, adj_mat, visited, group, threshold=threshold)
                    connected_groups_idx.append(group)
            return connected_groups_idx


        def DFS_recurse(idx, adj_mat, visited, connected_list, threshold):
            """Private: used by _get_connected_groups_init for recursion"""
            visited[idx] = True
            connected_list.append(idx)
            for i in range(len(visited)):
                if i == idx:
                    continue
                elif adj_mat[i][idx] <= threshold and not visited[i]:
                    _, _, visited, connected_list, _ = DFS_recurse(i, adj_mat,visited, connected_list, threshold)
            return idx, adj_mat, visited, connected_list, threshold

        connected_groups_idx = get_connected_groups_init(peaks, adj_mat, threshold)
        accepted_groups_idx = []
        for group_idx in connected_groups_idx:
            flag = False
            for i, j in combinations(group_idx, r=2):
                if adj_mat[i][j] > threshold*3:
                    group_vals = [peaks[i] for i in group_idx]
                    group_matrix = Peak.get_adjacency_matrix(group_vals, seq_len=200)
                    split_group_idx = get_connected_groups_init(group_vals, group_matrix, threshold//2)
                    split_group = [[group_idx[i] for i in group] for group in split_group_idx]
                    accepted_groups_idx.extend(split_group)
                    flag = True
                    break
            if not flag:
                accepted_groups_idx.append(group_idx)
            else:
                flag = False
                
        connected_groups_vals = [ [peaks[i] for i in idx_group] for idx_group in accepted_groups_idx ]
        return connected_groups_vals


    @staticmethod
    def get_middle(point_a: Union[int, "Peak"] , point_b: Union[int, "Peak"], curve_size: int = None) -> int:
        """
        Calculate the distance between Peak a and b on circular DNA.
        Returns a new Peak.middle (int) in the middle of a and b
        """
        a = point_a if not isinstance(point_a, Peak) else point_a.middle
        b = point_b if not isinstance(point_b, Peak) else point_b.middle
        seq_len = point_a.seq_len if curve_size is None and isinstance(point_a, Peak) else curve_size

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

        # def intersecting_windows(self, other: "Peak") -> bool:
        #     '''T|F wether self's window intersects with other's window'''
        #     if not isinstance(other, Peak):
        #         raise ValueError(f'other is a {type(other)}, must be a Peak object')
        #     # f = five_side, t = three_side
        #     f_t = Peak.calc_dist(self.five_side, other.three_side, self.seq_len)
        #     t_f = Peak.calc_dist(self.three_side, other.five_side, self.seq_len)
        #     f_f = Peak.calc_dist(self.five_side, other.five_side, self.seq_len)
        #     t_t = Peak.calc_dist(self.three_side, other.three_side, self.seq_len)

        #     a = f_t <= self.window_size // 2 or t_f <= self.window_size // 2
        #     b = f_f <= self.window_size // 2 and t_t <= self.window_size // 2
        #     return a or b

        if not isinstance(other, Peak):
            raise ValueError(f'other is a {type(other)}, must be a Peak object')
        x = Peak.calc_dist(self.middle, other.middle, self.seq_len)
        y = self.window_size//2 + other.window_size//2
        return x < y


    def contains_point(self, point: Union[int, "Peak"]) -> bool:
        '''T|F wether a point is within a Peak's window'''
        p = point.middle if isinstance(point, Peak) else point
        return self.window_size // 2 >= Peak.calc_dist(self.middle, p, self.seq_len)


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

    def __repr__(self):
        return f'Peak(middle={self.middle}, window_size={self.window_size})'

    def __str__(self):
        return str(self.middle)

    def __lt__(self, other: "Peak") -> bool:
        return self.middle < other.middle

    def __gt__(self, other: "Peak") -> bool:
        return self.middle > other.middle

    def __add__(self, other: "Peak") -> Union["Peak", int, float]:
        if isinstance(other, Peak):
            return Peak(self.middle + other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle + other

    def __sub__(self, other: "Peak") -> Union["Peak", int, float]:
        if isinstance(other, Peak):
            return Peak(self.middle - other.middle, self.seq_len, self.window_size)
        elif isinstance(other, int) or isinstance(other, float):
            return self.middle - other

    def __radd__(self, other: "Peak") -> Union["Peak", int, float]:
        return self.__add__(other)

    def __rsub__(self, other: "Peak") -> Union["Peak", int, float]:
        return self.__sub__(other)