"""Module for handing and processing disparity curves"""

from itertools import combinations

import numpy as np
from scipy.signal import find_peaks

from orca.Peak import Peak


def process_curve(curve: np.ndarray, mode: str, window_size: int) -> list[Peak]:
    '''Runs the given 1D-array (curve) through all processing functions for oriC identification. Returns its peaks'''
    init_peaks = [Peak(peak, curve.shape[0], window_size) for peak in detect_peaks(curve)]
    accepted_peaks = filter_peaks(curve, mode, init_peaks)
    peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

    single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
    merged_peaks = [Peak.from_calc_middle(to_merge[0], to_merge[1], curve.shape[0],
                                          window_size) for to_merge in peaks_to_merge]
    return single_peaks + merged_peaks


def curve_combinations(peaks_list: tuple[list["Peak"]], seq_len: int) -> list:
    '''Get every matched_peaks combination for x, y, and gc.'''
    oriC_locations_list = []
    for peaks_i, peaks_j in combinations(peaks_list, 2):
        matched_peaks = Peak.match_peaks(peaks_i, peaks_j)
        oriC_locations_list.append(
            [Peak.from_calc_middle(matches[0], matches[1], seq_len, peaks_list[0][0].window_size)
             for matches in matched_peaks]
        )
    return oriC_locations_list


def detect_peaks(curve: np.ndarray) -> np.ndarray:
    '''Calculates peaks of 1D-np.array and returns its indices.'''
    maxima, _ = find_peaks(curve, distance=curve.shape[0]//12)
    maxima = np.append(maxima, curve.argmax())
    minima, _ = find_peaks(np.negative(curve), distance=curve.shape[0]//12)
    minima = np.append(minima, curve.argmin())
    return np.unique(np.concatenate((maxima, minima), axis=0))


def filter_peaks(curve: np.ndarray, mode: str, peaks: list[Peak]) -> list[Peak]:
    '''
    Filters the given peaks based on the type of extreme it should be and the area around the peak in two steps.
    - Filter 1: Check if peaks are actually the extreme in their windows. This accounts for local extremes on the sides of slopes.
    - Filter 2: Check if any windows intersect the window of another peak (circular DNA) and reject the one that was less extreme (min or max).
    ### Input:
    - `curve`          : 1D-np.array
    - `peaks`          : list of Peaks of curve
    - `mode`           : 'max'|'min'. Which type of extreme do you want to find?
    ### Return:
    - `accepted_peaks` : list of peaks that passed both filters
    '''
    rejected_peaks = []
    rejected_peaks.extend(_filter_within_windows(curve, mode, peaks))
    rejected_peaks.extend(_filter_intersecting_windows(curve, mode, peaks))

    # Return list of peaks that passed both filters
    return list(set(peaks).difference(rejected_peaks))


def _filter_within_windows(curve: np.ndarray, mode: str, peaks: list[Peak]):
    """Filter 1: Check if peaks are actually the extreme in their windows"""
    rejected = []
    for peak in peaks:
        if peak.split:
            # Determine which side of 0 the peak.middle (m) is:
            #     _______0                             0_______
            #     5'   m    3'                      5'   m    3'
            # -----------|-----------  vs.  -----------|-----------
            comparator_win = curve[peak.five_side:] if (peak.middle > peak.five_side) else curve[:peak.three_side]
        else:
            comparator_win = curve[peak.five_side:peak.three_side]

        if mode == 'max' and comparator_win.max() > curve[peak.middle]:
            rejected.append(peak)
        elif mode == 'min' and comparator_win.min() < curve[peak.middle]:
            rejected.append(peak)
    return rejected


def _filter_intersecting_windows(curve, mode, peaks):
    rejected = []
    for peak_i, peak_j in combinations(peaks, 2):
        # Filter 1: Check if any windows intersecting the window of peak i
        if peak_i.intersecting_windows(peak_j):
            # Either peaks at the beginning and end of the DNA sequence or a peak found by sp.find_peaks() that is very close to the global min/max
            if mode == 'max' and curve[peak_i.middle] > curve[peak_j.middle]:
                rejected.append(peak_j)
            elif mode == 'min' and curve[peak_i.middle] < curve[peak_j.middle]:
                rejected.append(peak_j)
            else:
                rejected.append(peak_i)
    return rejected
