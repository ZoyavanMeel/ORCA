import os, csv, pickle
from itertools import combinations
from typing import TextIO, Union
from urllib.error import HTTPError, URLError

import numpy as np
from scipy.signal import find_peaks
from Bio import SeqIO, Entrez

from Peak import Peak


class FileHandler:
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
            raise ValueError(f'Unable to fetch accession: \'{accession}\' using \'{email}\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.') from BadRequest
        except URLError as NoConnection:
            raise ConnectionError('You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.') from NoConnection


    @staticmethod
    def parse_SeqRecord(record: SeqIO.SeqRecord, genes_of_interest: list[str]) -> dict:
        """Extracts the sequence and positions of the genes of interest from a SeqRecord."""
        seq_dict = dict()
        accession, version = tuple(record.id.split('.'))
        seq_dict.update({'accession': accession})
        seq_dict.update({'version': int(version)})
        seq_dict.update({'seq': record.seq})
        seq_dict.update({'seq_len': len(record.seq)})
        seq_dict.update({'gene_locations': []})
        seq_dict.update({'NCBI_oriC': []})

        for feature in record.features:
            # is this feature a coding sequence and a gene and is its name something we are looking for?
            if feature.type == 'CDS' and 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in genes_of_interest:
                gene_loc = Peak.from_edges(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
                seq_dict['gene_locations'].append( (feature.qualifiers['gene'][0], gene_loc) )
            # just in case this SeqRecord has an annotated oriC!
            if feature.type == 'rep_origin':
                oriC = Peak.from_edges(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
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



class CurveHandler:
    """Static class for handing and processing disparity curves"""
    @staticmethod
    def process_curve(curve: np.ndarray, mode: str, window_size: int) -> list:
        '''Runs the given 1D-array (curve) through all processing functions for oriC identification. Returns its peaks'''
        init_peaks = [Peak(peak, curve.shape[0], window_size) for peak in CurveHandler.detect_peaks(curve)]
        accepted_peaks = CurveHandler.filter_peaks(curve, mode, init_peaks)
        peaks_to_merge = Peak.get_peaks_to_merge(accepted_peaks)

        single_peaks = [x for x in accepted_peaks if not any(x in y for y in peaks_to_merge)]
        # merged_peaks = [Peak(to_merge[0].get_middle(to_merge[1]), curve.shape[0], window_size) for to_merge in peaks_to_merge]
        merged_peaks = [Peak.from_edges(to_merge[0], to_merge[1], curve.shape[0], window_size) for to_merge in peaks_to_merge]
        return single_peaks + merged_peaks
    

    @staticmethod
    def curve_combinations(peaks_list: tuple[list["Peak"]], seq_len: int) -> list:
        '''Get every matched_peaks combination for x, y, and gc.'''
        oriC_locations_list = []
        for peaks_i, peaks_j in combinations(peaks_list, 2):
            matched_peaks  = Peak.match_peaks(peaks_i, peaks_j)
            oriC_locations_list.append(
                # Peak( Peak.get_middle(matches[0], matches[1]), seq_len, peaks_list[0][0].window_size
                [Peak.from_edges(matches[0], matches[1], seq_len, peaks_list[0][0].window_size ) for matches in matched_peaks]
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
                if mode == 'max':
                    reject_middle = np.where( curve == min(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()
                elif mode == 'min':
                    reject_middle = np.where( curve == max(curve[peak_i.middle], curve[peak_j.middle]) )[0].tolist()

                if peak_i.middle in reject_middle:
                    rejected_peaks.append(peak_i)
                else:
                    rejected_peaks.append(peak_j)

        for peak in peaks:
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
