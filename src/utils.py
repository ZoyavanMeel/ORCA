import os, csv, pickle, socket
from itertools import combinations
from typing import TextIO, Union
from urllib.error import HTTPError, URLError
import matplotlib.pyplot as plt

import numpy as np
from scipy.signal import find_peaks
from Bio import SeqIO, Entrez

from Peak import Peak


class FileHandler:
    """Static class for handing biodata files and whatnot."""
    @staticmethod
    def fetch_file(accession: str, email: str, api_key: Union[str, None], rettype: str) -> TextIO:
        """Downloads the given file_type of the given accession in temporary memory"""
        Entrez.email = email
        if api_key is not None: Entrez.api_key = api_key
        try:
            return Entrez.efetch(db="nuccore", id=accession, rettype=rettype, retmode="text")
        except HTTPError as BadRequest:
            if api_key is not None:
                api_message_1 = f' with API_key \'{api_key}\''
                api_message_2 = ' and if your API_key is correctly typed and linked to the given email'
            else :
                api_message_1, api_message_2 = '', ''
            BadRequest.add_note(f'Unable to fetch accession: \'{accession}\' using \'{email}\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.')
            raise
        except URLError as NoConnection:
            NoConnection.add_note('You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.')
            raise


    @staticmethod
    def parse_SeqRecord(record: SeqIO.SeqRecord, genes_of_interest: list[str]) -> dict:
        """Extracts the sequence and positions of the genes of interest from a SeqRecord."""
        seq_dict = dict()
        accession, version = tuple(record.id.split('.'))
        seq_dict.update({'accession': accession})
        seq_dict.update({'version': int(version)})
        seq_dict.update({'seq': str(record.seq)}) # I only use the sequence to loop over once. strings are much faster for this.
        seq_dict.update({'seq_len': len(record.seq)})
        seq_dict.update({'gene_locations': []})
        seq_dict.update({'NCBI_oriC': []})

        for feature in record.features:
            # is this feature a coding sequence and a gene and is its name something we are looking for?
            if feature.type == 'CDS' and 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in genes_of_interest:
                gene_loc = Peak.from_calc_middle(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
                seq_dict['gene_locations'].append( (feature.qualifiers['gene'][0], gene_loc) )
            # just in case this SeqRecord has an annotated oriC!
            if feature.type == 'rep_origin':
                oriC = Peak.from_calc_middle(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
                seq_dict['NCBI_oriC'].append( ('oriC', oriC) )
        return seq_dict


    @staticmethod
    def save_gbk(accession: str, email: str, output_folder: str, api_key: str = None):
        '''
        Download the GenBank file of a given accession and save it into the output_folder.
        Name of the file will be: `{accession}_{version}.gbk`.
        Raises a `FileExistsError` if a file with the generated name already exists in the output folder.
        If version is not provided in the accession, then the function downloads the latest version.
        '''
        with FileHandler.fetch_file(accession, email, api_key, rettype="gbwithparts") as fh:
            # Quick search for the version of the accession that was downloaded.
            _acc, version = FileHandler.get_accession_from_gbk(fh)
            
            # Check if a file with the same name already exists
            if any(file.startswith(_acc + '_' + version + '.gbk') for file in os.listdir(output_folder)):
                fh.close()
                raise FileExistsError(f'\'{_acc}_{version}.gbk\' already exists in: {output_folder}')
            
            # Save contents to path
            file_path = os.path.join(output_folder, _acc + '_' + version + '.gbk')
            with open(file_path, 'w') as oh:
                oh.write(fh.read())
                oh.close()
            fh.close()


    @staticmethod
    def save_pkl(accession: str, email: str, output_folder: str, api_key: str = None):
        '''
        Download the GenBank file of a given accession, parses it with Biopython into a SeqRecord, and save it into the output_folder.
        Name of the file will be: `{accession}_{version}.pkl`.
        Raises a `FileExistsError` if a file with the generated name already exists in the output folder.
        If version is not provided in the accession, then the function downloads the latest version.
        '''
        with FileHandler.fetch_file(accession, email, api_key, rettype="gbwithparts") as fh:
            # Quick search for the version of the accession that was downloaded.
            _acc, version = FileHandler.get_accession_from_gbk(fh)
            
            # Check if a file with the same name already exists
            if any(file.startswith(_acc + '_' + version + '.pkl') for file in os.listdir(output_folder)):
                fh.close()
                raise FileExistsError(f'\'{_acc}_{version}.pkl\' already exists in: {output_folder}')
            
            # Parse gbk file
            seq_rec = SeqIO.read(fh, 'gb')
            
            # Save contents to path
            file_path = os.path.join(output_folder, _acc + '_' + version + '.pkl')
            with open(file_path, 'w') as oh:
                pickle.dump(seq_rec, oh)
                oh.close()
            fh.close()


    @staticmethod
    def get_accession_from_gbk(fh: TextIO) -> tuple[str, str]:
        """Reads the accession and version number from an open GenBank file."""
        acc_v = "."
        for line in fh:
            # Always appears around the top of the file.
            if "VERSION" in line:
                acc_v = [x for x in line.strip("\n").split(" ")][-1]
                break
        # Reset read to top of file.
        fh.seek(0)
        return tuple(acc_v.split("."))


    @staticmethod
    def merge_csvs(file_folder: str, merged_csv: str, fieldnames: list[str], length: int = -1, headers: bool = False):
        '''
        Merge multiple csvs into one csv.
        Arguments:
        - file_folder : path to folder with csvs that have to be merged
        - merged_csv  : name of the merged csv
        - fieldnames  : list of fieldnames for the merged_csv
        - length      : amount of rows in each single csv. -1, if the length of each single
                        csv is not known or not the same.
        - headers     : if the single csvs have headers or not
        '''
        file_list = os.listdir(file_folder)
        with open(merged_csv, 'w', newline='') as fh_out:
            writer = csv.writer(fh_out)
            writer.writerow(fieldnames)
            for file in file_list:
                with open(file_folder + file, 'r') as fh_in:
                    reader = csv.reader(fh_in)
                    for i, row in enumerate(reader):
                        if headers and i == 0:
                            pass
                        elif i == length:
                            break
                        else:
                            writer.writerow(row)



class Plotter:
    """Static class for plotting disparity curves."""
    @staticmethod
    def plot_Z_curve_3D(Z_curve: tuple[np.ndarray, np.ndarray, np.ndarray], name: str = None) -> None:
        """
        3D-plot function with name as title
        
        -----------------------------------------------------------------------
        Example
        >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
        >>> orca.find_oriCs()
        >>> plot_Z_curve_3D(Z_curve=[orca.x, orca.y, orca.z], name=orca.accession)
        """
        fig = plt.figure(figsize=(7,7))
        ax = plt.axes(projection='3d')
        ax.set_xlabel('Purine vs. Pyrimidine', fontsize=10, labelpad=15)
        ax.set_ylabel('Amino vs. Keto', fontsize=10, labelpad=15)
        ax.set_zlabel('Weak vs. Strong H-bond', fontsize=10, labelpad=15)

        x, y, z = Z_curve
        ax.plot3D(x, y, z, c='b', linewidth=0.8)
        if name is not None:
            ax.set_title(f'Z-Curve: {name}', fontsize=10, loc='center', pad=20)
        plt.show()


    @staticmethod
    def plot_x_curves(curves: tuple[np.ndarray], peaks: list[int], labels: list[str], name: str = None) -> None:
        """
        Plots up to 4 different y-axes onto a single figure. Ideal for displaying multiple disparity curves in a single plot.
        - `curves` : list of lists with y-axis values.
        - `peaks`  : list with indeces to plot onto the `curves`.
        - `labels` : list of names of each curve in `curves`.
        - `name`   : used in plot title.

        -----------------------------------------------------------------------
        Example
        >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
        >>> orca.find_oriCs()
        >>> plot_Z_curve_2D(curves=[orca.x, orca.y, orca.gc], peaks=orca.oriC_middles, labels=['$x_n$', '$y_n$', '$g_n$'], name=orca.accession)
        """
        x_len = str(len(curves[0]))
        if int(x_len[1]) <= 4:
            x_max = x_len[0] + str(int(x_len[1])+(5-int(x_len[1]))) + '0'*len(x_len[2:])
        else:
            x_max = str(int(x_len[0])+1) + '0'*len(x_len[1:])
        thing = int(x_max)//1000
        xthing = thing * 100
        ything = thing*2

        color_list = ['r', 'b', 'g', 'c']
        fig = plt.figure(figsize=(8,4))
        fig.subplots_adjust(right=0.75, bottom=0.25)
        base_ax = plt.axes()
        ax_list = [base_ax] + [base_ax.twinx() for i in range(len(curves) - 1)]

        offset = 1
        for axis in ax_list[1:]:
            axis.spines.right.set_position(("axes", offset))
            offset += 0.2

        good = False
        while not good:
            yticks_len = 0
            for i, ax in enumerate(curves):
                ubound, lbound = ything * round(max(ax)/ything), ything * round(min(ax)/ything)
                upper = ubound if max(ax) <= ubound else ubound + ything
                lower = lbound if min(ax) >= lbound else lbound - ything
                if len([x for x in range(lower, upper+ything, ything)]) > yticks_len:
                    yticks_len = len([x for x in range(lower, upper+ything, ything)])
            if yticks_len < 6:
                yticks_len *= 2
                ything = ything // 2
            else:
                good = True
                break

        handle_list = []
        for i, ax in enumerate(curves):
            peaks_y = np.asarray([ax[j] for j in peaks]) # y refers to the y-axis coordinates, not the y-curve
            ax_list[i].plot(range(len(ax)), ax, color=color_list[i], zorder=2, label=labels[i])
            ax_list[i].scatter(peaks, peaks_y, marker='o', c='k', zorder=3, label='$\it{oriC}$')
            ax_list[i].tick_params(axis ='y', colors=color_list[i])
            ax_list[i].ticklabel_format(axis='y', style='sci', useMathText=True)

            lbound = ything * round(min(ax)/ything)
            lower = lbound if min(ax) >= lbound else lbound - ything
            yticks = [lower + ything*j for j in range(yticks_len)]
            ax_list[i].set_yticks(yticks)
            ax_list[i].set_ylim(min(yticks), max(yticks))

            h, _ = ax_list[i].get_legend_handles_labels()
            handle_list.extend(h[:-1])
            oriC_handle = h[-1]
        handle_list.append(oriC_handle)

        if name is not None:
            base_ax.set_title(f'2D Z-Curve: {name}', fontsize=10,loc='center', pad=20)

        ubound= xthing * round(len(ax)/xthing)
        upper = ubound if len(ax) <= ubound else ubound + xthing
        xticks = [x for x in range(0, upper+xthing, xthing)]
        base_ax.set_xticks(xticks)

        base_ax.ticklabel_format(axis='x', style='sci', scilimits=(3,3), useMathText=True)
        base_ax.set_xlabel('Sequence length (bp)')
        base_ax.set_xlim(min(xticks), max(xticks))
        base_ax.grid(True, which='major', zorder=1)

        plt.legend(
            handles=handle_list,
            labels=labels + ['Prediction'],
            bbox_to_anchor=(0.12, -0.35, 0.75, .102),
            loc='center',
            ncol=len(curves)+1,
            mode="expand",
            borderaxespad=0.
        )
        plt.show()


    @staticmethod
    def plot_skew(skewArray: np.ndarray, peaks: list[int], name: str) -> None:
        """Plots single skew diagram and its peaks"""

        fig = plt.figure()
        ax1 = plt.axes()

        peaks_y = np.asarray([skewArray[i] for i in peaks])

        ax1.set_title(f'GC-skew: {name}', fontsize=10,loc='center', pad=20)
        ax1.plot(range(len(skewArray)), skewArray, 'r', zorder=1)
        ax1.scatter(peaks, peaks_y, marker='X', c='k', zorder=2)
        plt.show()



class CurveHandler:
    """Static class for handing and processing disparity curves"""
    @staticmethod
    def process_curve(curve: np.ndarray, mode: str, window_size: int) -> list[Peak]:
        '''Runs the given 1D-array (curve) through all processing functions for oriC identification. Returns its peaks'''
        init_peaks = [Peak(peak, curve.shape[0], window_size) for peak in CurveHandler.detect_peaks(curve)]
        accepted_peaks = CurveHandler.filter_peaks(curve, mode, init_peaks)
        peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

        single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
        merged_peaks = [Peak.from_calc_middle(to_merge[0], to_merge[1], curve.shape[0], window_size) for to_merge in peaks_to_merge]
        return single_peaks + merged_peaks
    

    @staticmethod
    def curve_combinations(peaks_list: tuple[list["Peak"]], seq_len: int) -> list:
        '''Get every matched_peaks combination for x, y, and gc.'''
        oriC_locations_list = []
        for peaks_i, peaks_j in combinations(peaks_list, 2):
            matched_peaks  = Peak.match_peaks(peaks_i, peaks_j)
            oriC_locations_list.append(
                [Peak.from_calc_middle(matches[0], matches[1], seq_len, peaks_list[0][0].window_size ) for matches in matched_peaks]
            )
        return oriC_locations_list


    @staticmethod
    def detect_peaks(curve: np.ndarray) -> np.ndarray:
        '''Calculates peaks of 1D-np.array and returns its indeces.'''
        maxima, _ = find_peaks( curve, distance=curve.shape[0]//12)
        maxima    = np.append(maxima, curve.argmax())
        minima, _ = find_peaks( np.negative(curve), distance=curve.shape[0]//12)
        minima    = np.append(minima, curve.argmin())
        return np.unique(np.concatenate( (maxima, minima), axis=0))
    

    @staticmethod
    def filter_peaks(curve: np.ndarray, mode: str, peaks: list[Peak]) -> list[Peak]:
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
                if mode == 'max' and curve[peak_i.middle] > curve[peak_j.middle]:
                    rejected_peaks.append(peak_j)
                elif mode == 'min' and curve[peak_i.middle] < curve[peak_j.middle]:
                    rejected_peaks.append(peak_j)
                else:
                    rejected_peaks.append(peak_i)

        for peak in peaks:
            # TODO: This second filter sucks wtf. Fix this...
            # Filter 2: Check if peaks are actually the extreme in their windows
            if peak.split:
                a, b = (peak.five_side, len(curve)-1), (0, peak.three_side)
                comparator_win = a if (peak.middle >= a[0] and peak.middle <= a[1]) else b
            else:
                comparator_win = (peak.five_side, peak.three_side)

            if mode == 'max' and np.max( curve[comparator_win[0]:comparator_win[1]] ) > curve[peak.middle]:
                    rejected_peaks.append(peak)
            elif mode == 'min' and np.min( curve[comparator_win[0]:comparator_win[1]] ) < curve[peak.middle]:
                    rejected_peaks.append(peak)

        # Create list of peaks that passed both filters
        rejected_peaks = set(rejected_peaks)
        accepted_peaks = [x for x in peaks if x not in rejected_peaks]
        return accepted_peaks
