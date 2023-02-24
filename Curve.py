import numpy as np
import scipy.signal as sp
from itertools import combinations
from Peak import Peak

class Curve:
    '''Object for handling curves, i.e. x, y, z, gc components.'''
    def __init__(self, curve: np.ndarray, mode: str, peaks: list[str] = None):
        self.curve = curve
        self.mode = mode
        self.peaks = [] if peaks is None else peaks


    def process_array(self, window_size: int) -> list:
        '''Runs the given 1D-array (curve) through all processing functions for oriC identification. Returns its peaks'''
        init_peaks = [Peak(peak, len(self.curve), window_size) for peak in self.detect_peaks()]
        accepted_peaks = self.filter_peaks(init_peaks)
        peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

        single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
        merged_peaks = [Peak(to_merge[0].get_middle(to_merge[1]),len(self.curve), window_size) for to_merge in peaks_to_merge]
        return single_peaks + merged_peaks
    

    def detect_peaks(self) -> np.ndarray:
        '''Calculates peaks of 1D-np.array and returns its indeces.'''
        maxima, _ = sp.find_peaks( self.curve, distance=len(self.curve)//12)
        maxima    = np.append(maxima, self.curve.argmax())
        minima, _ = sp.find_peaks( np.negative(self.curve), distance=len(self.curve)//12)
        minima    = np.append(minima, self.curve.argmin())
        return np.unique(np.concatenate( (maxima, minima), axis=0))
    

    def filter_peaks(self, peaks: list[Peak]) -> list:
        '''
        Filters the given peaks based on the type of extreme it should be and the area around the peak in two steps.
        - Filter 1: Check if any windows intersect the window of another peak (circular DNA).
        - Filter 2: Check if peaks are actually the extreme in their windows.\n
        Input:
        - `curve`          : 1D-np.array
        - `peaks`          : list of Peaks of curve
        - `mode`           : 'max'|'min'. Which type of extreme do you want to find?\n
        Return:
        - `accepted_peaks` : peaks that passed both filters
        '''
        rejected_peaks = []
        for peak_i, peak_j in combinations(peaks, 2):
            # Filter 1: Check if any windows intersecting the window of peak i
            if peak_i.intersecting_windows(peak_j):
                # Either peaks at the beginning and end of the DNA sequence or a peak found by sp.find_peaks() that is very close to the global min/max 
                if self.mode == 'max': reject_middle = np.where( self.curve == min(self.curve[peak_i.middle], self.curve[peak_j.middle]) )[0].tolist()
                elif self.mode == 'min': reject_middle = np.where( self.curve == max(self.curve[peak_i.middle], self.curve[peak_j.middle]) )[0].tolist()
                rejected_peaks.append(peak_i) if peak_i.middle in reject_middle else rejected_peaks.append(peak_j)

        for peak in peaks:
            # Filter 2: Check if peaks are actually the extreme in their windows
            if peak.split:
                a, b = (peak.five_side, len(self.curve)-1), (0, peak.three_side)
                comparator_win = a if (peak.middle >= a[0] and peak.middle <= a[1]) else b
            else:
                comparator_win = (peak.five_side, peak.three_side)

            if self.mode == 'max' and np.max( self.curve[comparator_win[0]:comparator_win[1]] ) > self.curve[peak.middle]:
                    rejected_peaks.append(peak)
            elif self.mode == 'min' and np.min( self.curve[comparator_win[0]:comparator_win[1]] ) < self.curve[peak.middle]:
                    rejected_peaks.append(peak)

        # Create list of peaks that passed both filters
        rejected_peaks = set(rejected_peaks)
        accepted_peaks = [x for x in peaks if x not in rejected_peaks]
        return accepted_peaks