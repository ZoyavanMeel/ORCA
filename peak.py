from typing import Union, Tuple, List
from itertools import combinations

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


    @staticmethod
    def calc_dist( a: int, b: int, curve_size: int) -> int:
        '''Calculate the distance of self to other in bp'''
        dist_1 = max(a, b) - min(a, b)
        dist_2 = min(a, b) + curve_size-1 - max(a, b)
        return min(dist_1, dist_2)


    @staticmethod
    def get_peaks_to_merge(peaks: list) -> List[Tuple["Peak", "Peak"]]:
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
        '''T|F wether self's window intersects with other's window'''
        if not isinstance(other, Peak):
            raise ValueError(f'other is a {type(other)}, must be a Peak object')
        # f = five_side, t = three_side
        f_t, t_f = Peak.calc_dist(self.five_side, other.three_side, self.seq_len), Peak.calc_dist(self.three_side, other.five_side, self.seq_len)
        f_f, t_t = Peak.calc_dist(self.five_side, other.five_side, self.seq_len), Peak.calc_dist(self.three_side, other.three_side, self.seq_len)
        a, b = f_t <= self.window_size // 2 or t_f <= self.window_size // 2, f_f <= self.window_size // 2 and t_t <= self.window_size // 2
        return a or b


    def contains_point(self, point: Union[int, "Peak"]) -> bool:
        '''T|F wether a point is within a Peak's window'''
        p = point.middle if isinstance(point, Peak) else point
        return self.window_size // 2 >= Peak.calc_dist(self.middle, p, self.seq_len)


    # Undefined __eq__ to remain hashable; only added the dunders I needed
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