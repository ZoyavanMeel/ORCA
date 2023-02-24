# Libraries
import warnings
import numpy as np

from itertools import combinations, product
from typing import Tuple, Generator

# Self-made modules
from Peak import Peak
from SequenceHandler import SequenceHandler
from GeneHandler import GeneHandler
import plotter_functions as pf


class ORCA:
    def __init__(self, **kwargs):
        self.args = ORCA._check_variables(**kwargs)


    @staticmethod
    def _check_variables(**kwargs):
        """Handle the user input parameters and throw errors and warnings if they are not valid."""
        defaults = {
            'accession'         : None,
            'email'             : None,
            'api_key'           : None,
            'genome_fasta'      : None,
            'genes_fasta'       : None,
            'dnaa_boxes'        : ['TTATACACA', 'TTATTCACA', 'TTATCCACA', 'TTATGCACA'],
            'genes_of_interest' : ['dnaA', 'dnaN'],
            'windows'           : [0.01, 0.03, 0.05],
            'max_group_spread'  : 0.05,
            'max_mismatches'    : 0,
            'model'             : None,
            'show_info'         : False,
            'show_plot'         : False
        }

        # Update defaults with provided values
        args = {**defaults, **kwargs}

        # Fasta retrieval checks
        if args['accession'] is not None and not isinstance(args['accession'], str):
            raise TypeError('accession must be a string.')
        if args['email'] is not None and not isinstance(args['email'], str):
            raise TypeError('email must be a string.')
        if args['api_key'] is not None and not isinstance(args['api_key'], str):
            raise TypeError('api_key must be a string.')

        # File format checks
        if isinstance(args['genome_fasta'], str) and not (args['genome_fasta'][-5:] == 'fasta' or args['genome_fasta'][-3:] == 'fna'):
            raise ValueError('\'genome_fasta\' does not have extension \'.fasta\' or \'.fna\'. Must be the path to a FASTA-file.')
        if isinstance(args['genes_fasta'], str) and not (args['genes_fasta'][-5:] == 'fasta' or args['genes_fasta'][-3:] == 'fna'):
            raise ValueError('\'genes_fasta\' does not have extension \'.fasta\' or \'.fna\'. Must be FASTA-file.')
        
        # Variable combination checks
        if args['genome_fasta'] is None and args['genes_fasta'] is None and args['accession'] is None:
            raise ValueError('Did not provide files to read or accession to fetch.')
        if args['genome_fasta'] is not None and args['accession'] is not None:
            warnings.warn('Provided both a fasta to read and an accession to fetch. Will ignore given accession and use accession from \'genome_fasta\'.')
        if args['accession'] is not None and args['email'] is None:
            raise ValueError('Did not provide a email adress for fetching the accession.\n\tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/\n\tCreate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/')
        if args['genome_fasta'] is not None and args['genes_fasta'] is None and args['email']  is None:
            raise ValueError('Only provided \'genome_fasta\'. Will have to fetch \'genes_fasta\', but you provided no \'email\'.')
        return args


    @staticmethod
    def _merge_oriCs(curve_size: int, groups: list, window_size: int) -> Tuple[list, list]:
        '''Finds the average index of a group and returns those values. groups is a nested-list'''
        mutable = sorted( groups, key=lambda x:len(x), reverse=True )
        total_pot_oriCs = len( [y for x in mutable for y in x] )
        oriCs, Z_scores = [], []

        group: list[Peak]
        for group in mutable:
            group.sort()
            for i in range(len(group)):
                if (group[-1] - group[i]).middle >= (curve_size-1)/2:
                    group[i].middle += curve_size-1
            avg_val = sum(group)//len(group)
            if avg_val > curve_size-1:
                avg_val -= curve_size-1
            oriCs.append(Peak(avg_val, curve_size, window_size))
            Z_scores.append(len(group)/total_pot_oriCs)
        return oriCs, Z_scores


    @staticmethod
    def get_adj_mat(peaks_a: list[Peak], peaks_b: list[Peak] = None, seq_len: int = None) -> np.ndarray:
        """
        Gets adjacency matrix for given list of `Peak` or `int` objects.
        The matrix can be between a list and the elements of itself or between the elements of two lists.
        All elements in each list must be of the same type.
        Elements in `peaks_a` must be the same type as those in `peaks_b`.
        If elements are integers, `seq_len` must be provided.
        """

        def get_adj_mat_setup(peaks_a: list, peaks_b: list = None, seq_len: int = None) -> tuple[bool, np.ndarray, Generator]:
            '''Check input variables for get_adj_mat and initialise adj_mat and iterator'''
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
            return are_integers, adj_mat, iterator

        # Error handling and variable initialisation
        are_integers, adj_mat, iterator = get_adj_mat_setup(peaks_a, peaks_b, seq_len)

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
                    group_matrix = ORCA.get_adj_mat(group_vals, seq_len=200)
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


    def find_oriCs(self):
        '''
        Locates potential oriCs on circular bacterial chromosomes based on Z-curve and GC-skew analysis, dnaA box analysis, and dnaA/dnaN gene locations.
        Three default window_sizes are used: 1, 3 and 5 % of the total genome length. See the README-file in the [GitHub repository](https://github.com/ZoyavanMeel/ORCA/)
        for more information or consult ORCA.pdf for more extensive results and analyses of ORCA.

        This function either reads a given FASTA and genes_fasta or fetches them using a given accession directly from the NCBI database and calculates its oriC.

        **kwargs:
        - `accession`           : Accession number of the sequence to fetch. If no version is provided, will fetch the newest version.
        - `email`               : Email adress of your NCBI account
        - `api_key`             : API Key for downloading from the NCBI database (E-Utils).
                                Optional (as of 2022/07/10), but necessary if fetching at 10 requests per second or more.
        - `genome_fasta`        : FASTA-file with circular bacterial DNA
        - `genes_fasta`         : FASTA-file with gene info in the same format as when acquired using `E-Utils(db='nuccore', rettype='fasta_cds_na')`
        - `dnaa_boxes`          : If None, will use the [consensus DnaA-box](https://doi.org/10.1093/bib/bbn031): `TTAT(A|T|C|G)CACA`.
                                Else, provide a list of 9 base strings. See the `get_dnaa_boxes` function in `helper_functions.py` for some more examples of dnaA-boxes.
                                Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
        - `max_mismatches`      : Maximum allowed mismatches before a 9-mer is considered to fit the dnaa_box. Recommended: 0; recommended max: 2.
        - `genes_of_interest`   : List of gene names to look for in `genes_fasta`.
        - `max_group_spread`    : Maximum spread a group can have when looking for connected groups.
        - `show_info`           : If True, prints info of ALL found oriCs. Good and bad.
        - `show_plot`           : If True, shows plot of ALL found oriCs. Good and bad. Should not be used for analysis -> Make a separate plot for the best oriCs.

        Return:
        - `properties`          : Dictionary with properties of all oriC-like regions.
                                NOTE: oriCs are NOT sorted by importance. Recommended way to rank: learning machine decision.
            - `'oriCs`          : List of Peak-objects. Each Peak has a index position on the given sequence and scores based on the analyses (see: ORCA.pdf).
            - `'dnaA_boxes'`    : Dictionary with dnaA-box 9-mers as keys and lists of position indices on the given DNA as values.
                                The indices refer to the position of the 5th base in the 9-mer.
            - etc.
        '''
        seq = SequenceHandler(self.args)

        # Finding potential oriCs based on Z-curve
        peaks = []
        for fraction in self.args['windows']:
            window_size = int(seq.length * fraction)
            peaks_x  = seq.x.process_array(window_size=window_size)
            peaks_y  = seq.y.process_array(window_size=window_size)
            peaks_gc = seq.gc.process_array(window_size=window_size)
            peaks.extend( [j for i in seq.curve_combinations( (peaks_x, peaks_y, peaks_gc) ) for j in i] )
        
        ## Finding connected components in undirected graph with a depth-first search to merge Z-curve oriCs
        matrix_pot_oriCs = ORCA.get_adj_mat(peaks)
        connected_groups = ORCA.get_connected_groups(peaks, matrix_pot_oriCs, int(seq.length*self.args['max_group_spread']))
        oriCs, Z_scores  = ORCA._merge_oriCs(seq.length, connected_groups, window_size=int(seq.length*self.args['windows'][-1]))