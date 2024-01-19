# Libraries
import os
import pickle
import re
import warnings
from itertools import combinations, product
from typing import Generator

import numpy as np
from Bio import Seq, SeqIO
from sklearn.base import is_classifier

# Self-made modules
import BioFile
import CurveProcessing
import Plotter
from Peak import Peak


class ORCA:
    """
    ORCA: Origin of Chromosomal Replication Assessment. Class for predicting the location of the origin of replication (ORC/oriC) of circular bacterial DNA.
    For further instruction, consult the README on [GitHub](https://github.com/ZoyavanMeel/ORCA/).

    Accepted **kwargs:
    ------------------
    - `dnaa_boxes`:
        - If None, will use the [consensus DnaA-box](https://doi.org/10.1093/bib/bbn031): `TTAT(A|T|C|G)CACA`.
        Else, provide a list of 9 base strings. See the `get_dnaa_boxes_with_mismatches` function for some more examples of dnaA-boxes.
        Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
    - `max_mismatches`:
        - Maximum allowed mismatches allowed in a dnaa_box for it still to be read as such. Recommended max: 2.
    - `genes_of_interest`:
        - List of gene names to consider as 'oriC-proximal' and use for helping estimate the location of the oriC. This parameter is case insensitive.
    - `max_point_spread`:
        - Maximum distance between points in a group can have when looking for connected groups across the disparity curves.
        Default is 5 % of the total chromosome length.
    - `windows`:
        - The windows around peaks of skew curves to consider. Defaults are 1, 3, and 5 % of the total chromosome length. ORCA checks each of the given windows.
    - `model`:
        - A fitted scikit-learn classifier. Recommended to use the one provided on [GitHub](https://github.com/ZoyavanMeel/ORCA/).
    """

    def __init__(self):
        """
        Empty constructor. Use: `from_gbk()`, `from_pkl()`, `from_accession()` or `from_string()` for instantiating ORCA objects.
        """
        pass

    def __make(self, **kwargs):
        """Private: used to check the validity of user-given keyword arguments and set them as instance variables."""
        user_args = ORCA.setup_variables(**kwargs)
        for key, value in user_args.items():
            setattr(self, key, value)

    @classmethod
    def from_string(cls, sequence: str, gene_locations: list[tuple[str, int, int]] = None, accession: str = None, **kwargs) -> "ORCA":
        """
        Minimal constructor: provide a string to analyse and whatever else you want.

        Parameters:
        ----------
        - `sequence`        : String representation of the DNA circular sequence to analyse.
        - `gene_locations`  : Optional, list of tuples in the following format: `(gene_name, start_position, end_position)`.
        - `accession`       : Optional, used for naming, else: 'Custom'.
        - `**kwargs`        : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        --------
        - `ORCA` : an ORCA object.

        --------------------------
        Example:
        --------

        >>> # Bio.Seq objects are also fine, they just take longer to loop over.
        >>> DNA = "ATCGATCGATCGATACGATGTGCTAGCTACTGATCGATCGACAGACTGCTAGCGATCCTCGA"
        >>> gene_locs = [("dnaA", 10, 30), ("dnaN", 40, 50)]
        >>> # To instantiate the ORCA object call:
        >>> orca = ORCA.from_string(DNA, gene_locs)
        >>> # To find oriCs:
        >>> orca.find_oriCs()
        """
        orca = cls()
        orca.__make(**kwargs)

        orca.seq = sequence
        orca.seq_len = len(sequence)
        orca.version = 0

        orca.accession = "Custom" if accession is None else accession
        orca.gene_locations = [] if gene_locations is None \
            else [(gene[0], Peak.from_calc_middle(gene[1], gene[2], len(sequence), 0)) for gene in gene_locations]

        return orca

    @classmethod
    def from_gbk(cls, path: str, **kwargs) -> "ORCA":
        """
        Instatiates an ORCA object from a GenBank (gbk) file.

        Parameters:
        ----------
        - `path`     : string of path to the GenBank file to analyse.
        - `**kwargs` : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        --------
        - `ORCA` : an ORCA object with properly loaded parameters.

        --------------------------
        Example:
        -------

        >>> # To instantiate the ORCA object call:
        >>> orca = ORCA.from_gbk("path/to/file.gbk")
        >>> # To find oriCs:
        >>> orca.find_oriCs()
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} does not exist.")

        orca = cls()
        orca.__make(**kwargs)

        file_args = BioFile.parse_SeqRecord(SeqIO.read(path, 'gb'), orca.genes_of_interest)
        for key, value in file_args.items():
            setattr(orca, key, value)
        return orca

    @classmethod
    def from_pkl(cls, path: str, **kwargs) -> "ORCA":
        """
        Instatiates an ORCA object from a pickled Biopython SeqRecord.

        Parameters:
        -----------
        - `path`     : string of path to the pickled SeqRecord of a GenBank file.
        - `**kwargs` : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        --------
        - `ORCA` : an ORCA object with properly loaded parameters.

        Raises:
        -------
        - `FileNotFoundError` : in case the file is not found at `path` 
        - `EOFError`          : in case of faulty/empty Pickle file.
        - `ValueError`        : in case the pickled object is not a SeqRecord.

        --------------------------
        Example:
        --------

        >>> # To instantiate the ORCA object call:
        >>> orca = ORCA.from_pkl("path/to/file.pkl")
        >>> # To find oriCs:
        >>> orca.find_oriCs()
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} does not exist.")

        orca = cls()
        orca.__make(**kwargs)

        with open(path, "rb") as fh:
            try:
                record = pickle.load(fh)
            except EOFError as eof:
                eof.add_note(
                    f"Something is wrong with your Pickle file. It is probably empty. ORCA was unable to load: \'{path}\'.")
                fh.close()
                raise
            fh.close()
            if not isinstance(record, SeqIO.SeqRecord):
                raise ValueError(f"Pickled object must be a SeqRecord object, but was a \'{type(record)}\' instead.")
            file_args = BioFile.parse_SeqRecord(record, orca.genes_of_interest)
            for key, value in file_args.items():
                setattr(orca, key, value)
            return orca

    @classmethod
    def from_accession(cls, accession: str, email: str, api_key: str = None, **kwargs) -> "ORCA":
        """
        Instatiates an ORCA object from a given accession number.

        Parameters:
        -----------
        - `accession` : accession to analyse. ORCA will fetch and parse the required file.
        - `email`     : your email address. Always tell NCBI who you are.
        - `api_key`   : Optional. Only necessary if you are doing more than 10 requests per second.
        - `**kwargs`  : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        --------
        - `ORCA` : an ORCA object with properly loaded parameters.

        --------------------------
        Example:
        --------

        >>> # To instantiate the ORCA object call:
        >>> orca = ORCA.from_accession("NC_example.3", "your@email.com")
        >>> # To find oriCs:
        >>> orca.find_oriCs()
        """
        orca = cls()
        orca.__make(**kwargs)

        if email is None:
            raise ValueError(f"""
                NCBI requires you to specify an email address to make use of their services.
                \tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/
                \tGenerate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/ (optional)
            """)
        with BioFile.fetch_file(accession, email, rettype='gbwithparts', api_key=api_key) as handle:
            file_args = BioFile.parse_SeqRecord(SeqIO.read(handle, 'gb'), orca.genes_of_interest)
            handle.close()
            for key, value in file_args.items():
                setattr(orca, key, value)
            return orca

    @staticmethod
    def setup_variables(**kwargs) -> dict:
        """
        Handle the user input parameters and throw errors and warnings if they are not valid. 
        Check the ORCA class documentation for valid keyword arguments.
        """

        defaults = {
            'dnaa_boxes': ['TTATACACA', 'TTATTCACA', 'TTATCCACA', 'TTATGCACA'],
            'genes_of_interest': ['dnaA', 'dnaN'],
            'windows': [0.01, 0.03, 0.05],
            'max_point_spread': 0.05,
            'max_mismatches': 0,
            'model': None,
        }

        # Update defaults with provided values
        args = {key: kwargs[key] if key in kwargs else defaults[key] for key in defaults}

        # Simple type enforcing. Overriding default parameters is okay, but provide empty lists if you do.
        assert args['model'] is None or is_classifier(args['model']), 'model must be a scikit-learn classifier.'
        assert isinstance(args['windows'], list), 'windows must be a list of integers.'
        assert isinstance(args['dnaa_boxes'], list), 'dnaa_boxes must be a list of strings.'
        assert isinstance(args['max_mismatches'], int), 'max_mismatches must be an integer.'
        assert isinstance(args['max_point_spread'], float), 'max_point_spread must be an integer.'
        assert isinstance(args['genes_of_interest'], list), 'genes_of_interest must be a list of strings.'

        # Setup all possible dnaA-box regions
        args['input_dnaa_boxes'] = args['dnaa_boxes'].copy()
        args['dnaa_boxes'] = ORCA.get_dnaa_boxes_with_mismatches(args['input_dnaa_boxes'], args['max_mismatches'])

        return args

    @staticmethod
    def get_dnaa_boxes_with_mismatches(box_list: list, max_mismatches: int = 2) -> list:
        """
        Standard parameters: Get all unique dnaa-box 9-mers and reverse complements that allow for 0, 1, or 2 mismatches.
        Sources as comments in function. Check source code for more suggested consensus boxes.

        Parameters:
        - `box_list`       : list with strings of dnaa-boxes to use. Make sure all given boxes are 9 bases long.
        - `max_mismatches` : maximum allowed mismatches in dnaa-box that is still seen as an accepted dnaa-box.
            E.g. 2 allows all kmers that have 0, 1 or 2 mismatches with the dnaa-box.\n
        Return:
            `dnaa_boxes`   : set of 9-mers matching the given parameters.

        ------------------------------------------------------------
        Useful acticles about DnaA(-boxes):
        - https://doi.org/10.1101/cshperspect.a012922
        - https://doi.org/10.1093/nar/gkr832
        - https://doi.org/10.3389/fmicb.2018.00319
        - https://doi.org/10.1046/j.1365-2958.1996.6481362.x.

        ------------------------------------------------------------

        |                                Sequence                                      |                    DOI                       |  Year  |                     Notes                         |
        | ---------------------------------------------------------------------------- | -------------------------------------------- | ------ | ------------------------------------------------- |
        | 'TTAT(A\|C)CA(A\|C)A'                                                        | https://doi.org/10.1016/0092-8674(84)90284-8 |  1984  | The first consensus                               |
        | 'TGTG(G\|T)ATAAC'                                                            | https://doi.org/10.1016/0022-2836(85)90299-2 |  1985  | Matsui-box                                        |
        | 'TTAT(A\|T\|C\|G)CACA'                                                       | https://doi.org/10.1093/dnares/dsm017        |  2007  | In (roughly) all eubacteria, not just B. subtilis |
        | 'TTATCCACA', 'TTTTCCACA', 'TTATCAACA', 'TCATTCACA', 'TTATACACA', 'TTATCCAAA' | https://doi.org/10.1093/emboj/16.21.6574     |  1997  | Affinity study                                    |
        | '(T\|C)(T\|C)(A\|T\|C)T(A\|C)C(A\|G)(A\|C\|T)(A\|C)'                         | https://doi.org/10.1007/BF00273584           |  1991  | Only in E. coli K12. Do not use.                  |
        """

        if max_mismatches == 0:
            return box_list

        def generate_mismatched_strings(string: str, mismatches: int = 2) -> Generator:
            string_len = len(string)
            letters = 'ACGT'

            for indices in combinations(range(string_len), mismatches):
                for replacements in product(letters, repeat=mismatches):
                    skip = False
                    for i, a in zip(indices, replacements):
                        if string[i] == a:
                            skip = True
                    if skip:
                        continue

                    keys = dict(zip(indices, replacements))
                    yield ''.join([string[i] if i not in indices else keys[i] for i in range(string_len)])

        # Get all dnaa-boxes as strings
        boxes = box_list.copy()
        for box in box_list:
            if re.search(r'[^ATCG]', box) is not None:
                raise ValueError(
                    f'\n\tInput string: \'{box}\' contains forbidden characters. Only use the four nucleotide bases: A, T, C, and G.')
            # Also check reverse complement for box on other strand
            boxes.append(str(Seq.Seq(box).reverse_complement()))

        # Get all unique strings while allowing for max. 2 mismatches.
        mismatch_boxes = list(set(boxes))
        for box in boxes:
            for i in range(abs(max_mismatches)):
                mismatch_boxes.extend(generate_mismatched_strings(box, i+1))
        return mismatch_boxes

    def calculate_disparity_curves(self) -> None:
        '''
        Z-curve and GC-skew calculation and k-mer indexing. In one function so only one iteration of the sequence is necessary.\n
        Sets:
            `x`, `y`, `z`, `gc` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
            `dnaa_dict`         : Dictionary of middle indexes of dnaa-boxes in `seq`
        '''
        k = len(self.dnaa_boxes[0])
        x, y, z, gc = [], [], [], []
        a, c, t, g = 0, 0, 0, 0

        raw_dict = {}

        for i in range(self.seq_len):
            base = self.seq[i]
            if base == "A":
                a += 1
            elif base == "C":
                c += 1
            elif base == "G":
                g += 1
            elif base == "T":
                t += 1

            gc.append(g - c)
            x.append((a + g) - (c + t))  # Purine vs Pyrimidine
            y.append((a + c) - (g + t))  # Amino vs Keto
            z.append((a + t) - (c + g))  # Weak vs Strong Hydrogen Bonds

            kmer = self.seq[i:i+k] if i <= self.seq_len - k else self.seq[i:] + self.seq[:k-(self.seq_len-i)]
            if kmer in self.dnaa_boxes:
                mid = i+k//2+1 if i+k//2+1 <= self.seq_len else i+k//2+1 - self.seq_len
                try:
                    raw_dict[kmer].append(mid)
                except KeyError:
                    raw_dict[kmer] = [mid]

        (self.x, self.y, self.z, self.gc, self.dnaa_dict) = (
            np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(gc), raw_dict)

    def analyse_disparity_curves(self) -> list[Peak]:
        '''
        Analyse the three disparity curves related to finding the oriC of the sequence.
        The Z-curve can only be analysed after calling `calc_disparities()`.
        '''
        peaks = []
        for fraction in self.windows:
            window_size = int(self.seq_len * fraction)
            peaks_x = CurveProcessing.process_curve(self.x, 'min', window_size=window_size)
            peaks_y = CurveProcessing.process_curve(self.y, 'max', window_size=window_size)
            peaks_gc = CurveProcessing.process_curve(self.gc, 'min', window_size=window_size)
            peaks.extend([j for i in CurveProcessing.curve_combinations(
                (peaks_x, peaks_y, peaks_gc), self.seq_len) for j in i])
        return peaks

    def calculate_Z_scores(self, peaks: list[Peak]) -> tuple[list[Peak], list[int]]:
        '''Finding connected components in undirected graph with a depth-first search to merge Z-curve oriCs'''

        adj_mat = Peak.get_adjacency_matrix(peaks)
        connected_groups = Peak.select_connected_groups(peaks, adj_mat, int(self.seq_len*self.max_point_spread))
        connected_groups.sort(key=lambda x: len(x), reverse=True)

        total_pot_oriCs = len([y for x in connected_groups for y in x])
        oriCs, Z_scores = [], []
        window_size = int(self.seq_len*self.windows[-1])

        for group in connected_groups:
            group.sort()
            for i in range(len(group)):
                if (group[-1] - group[i]).middle >= (self.seq_len-1)/2:
                    group[i].middle += self.seq_len-1
            avg_val = sum(group)//len(group)
            if avg_val > self.seq_len-1:
                avg_val -= self.seq_len-1
            oriCs.append(Peak(avg_val, self.seq_len, window_size))
            Z_scores.append(len(group)/total_pot_oriCs)
        return oriCs, Z_scores

    def calculate_D_scores(self, peaks: list[Peak]) -> list[Peak]:
        '''Process the location of dnaa_boxes and rank potential oriCs based on most surrounding dnaa_boxes.'''
        if len(self.dnaa_dict.keys()) == 0:
            return [0] * len(peaks)
        contains_boxes = []
        all_pos = [pos for pos_list in self.dnaa_dict.values() for pos in pos_list]
        for oriC in peaks:
            count = 0
            for pos in all_pos:
                if oriC.contains_point(pos):
                    count += 1
            contains_boxes.append(count)
        return [x/sum(contains_boxes) if sum(contains_boxes) != 0 else 0 for x in contains_boxes]

    def calculate_G_scores(self, peaks: list[Peak]) -> list[int]:
        '''Process the location of the genes of interest and rank the potential oriCs based on how close they are to these genes'''
        if len(self.gene_locations) == 0:
            return [0] * len(peaks)
        gene_peaks = [name_loc_tuple[1] for name_loc_tuple in self.gene_locations]
        matrix = Peak.get_adjacency_matrix(peaks, gene_peaks)
        if np.min(matrix) == np.max(matrix):
            return [1] * matrix.shape[1]
        norm_mat = (matrix - np.min(matrix)) / (np.max(matrix) - np.min(matrix))
        return [1 - x for x in np.mean(norm_mat, axis=1)]

    def __check_run_properly(self, attr: str, msg: str) -> None:
        """Used internally to check if the ORCA object has the proper attributes it needs at certain stages."""
        try:
            getattr(self, attr)
        except AttributeError as NotRunYet:
            NotRunYet.add_note(msg)
            raise

    def pretty_print_results(self) -> None:
        """Print the main ORCA results in a pretty way."""
        self.__check_run_properly("oriC_middles", "Cannot print results.")

        spacer_1 = len('oriC_middles')+1
        spacer_2 = len(str(max(self.oriC_middles)))+2

        print(f'{"Accession":<{spacer_1}}: {self.accession}')

        attrs = ['predictions', 'Z_scores', 'G_scores', 'D_scores', 'oriC_middles']
        for attr in attrs:
            if isinstance(getattr(self, attr)[0], float):
                results = ''.join(f'{round(x, spacer_2-4): >{spacer_2}}' for x in getattr(self, attr))
            else:
                results = ''.join(f'{str(x): >{spacer_2}}' for x in getattr(self, attr))
            print(f'{attr:<{spacer_1}}:{results}')
        print("The best-scoring potential oriC was found at:", self.oriC_middles[self.best_oriC_idx], "bp.")

    def plot_oriC_curves(self) -> None:
        """
        Plot curves relevant to ORCA's analysis as well as all found oriCs.
        To be called after calling `find_oriCs` on an ORCA object.
        See utils.Plotter for more plotting options.

        --------------------------------
        Example:

        >>> orca = ORCA.from_accession("NC_example.3", "your@email.com")
        >>> orca.find_oriCs()
        >>> orca.plot_oriC_curves()

        The same as calling:

        >>> orca = ORCA.from_accession('NC_000913', 'example@email.com')
        >>> orca.find_oriCs()
        >>> plot_curves(curves=[orca.x, orca.y, orca.gc], peaks=orca.oriC_middles, labels=['$x_n$', '$y_n$', '$GC_n$'], name=orca.accession)
        """
        self.__check_run_properly("oriC_middles", "Cannot plot curves.")
        Plotter.plot_curves([self.x, self.y, self.gc], ['$x_n$', '$y_n$', '$GC_n$'], self.oriC_middles, self.accession)

    def _set_best_pot_oriC_idx(self) -> None:
        """
        Set the best performing potential oriC based on the model's prediction as well as a tie-breaking method.

        Initially the potential oriC is chosen in which the model has the most "confidence", i.e. which one does
        the model think is the most likely to be the True origin. This is determined by the highest probability of
        "correctness" as calculated by the `predict_proba` function of the classifier.

        In the case of a tie, we look at the highest score of the most important feature (Z-score) between the
        tied potential oriCs. If this is also results in some form of tie, we check the next one, etc.
        """

        # set matrix columns to put them in order of feature importance: [prediction, Z-score, D-score, G-score]
        if len(self.predictions) == 0 or self.predictions[0] is None:
            mat = np.asarray([self.Z_scores, self.D_scores, self.G_scores]).T
        else:
            mat = np.asarray([self.predictions, self.Z_scores, self.D_scores, self.G_scores]).T
        max_idx_options = [j for j in range(mat.shape[0])]

        for i in range(mat.shape[1]):
            max_idx_in_curr_col = np.where(mat[max_idx_options, i] == np.max(mat[max_idx_options, i]))[0].tolist()
            max_idx_options = [max_idx_options[j] for j in max_idx_in_curr_col]
            if len(max_idx_options) == 1:
                break
        self.best_oriC_idx = max_idx_options[0]

    def find_oriCs(self, show_info: bool = False, show_plot: bool = False) -> None:
        '''
        Locates potential oriCs on circular bacterial chromosomes based on Z-curve and GC-skew analysis, dnaA box analysis, and dnaA/dnaN gene locations.
        Three default window_sizes are used: 1, 3 and 5 % of the total genome length. See the README-file in the [GitHub repository](https://github.com/ZoyavanMeel/ORCA/)
        for more information or consult ORCA.pdf for more extensive results and analyses of ORCA.

        Parameters:
        -----------
        - `show_info`:
            - If True, prints info of ALL found oriCs. Good and bad.

        - `show_plot`:
            - If True, shows plot of ALL found oriCs. Good and bad. Should not be used for analysis -> Make a separate plot for the best oriCs.

        Results
        -------
        Adds new instance variables to the ORCA object based on the analysis.

        All attributes after calling `find_oriCs`:
        - `accession`:
            - `str` Accesssion number.
        - `version`:
            - `int` Accession version number.
        - `input_dnaa_boxes`:
            - `list[str]` List of user-given dnaA-boxes to look for.
        - `dnaa_boxes`:
            - `list[str]` List of all possible accepted dnaA-boxes based on the given dnaA-boxes and the allowed `max_mismatches`.
        - `dnaa_dict`:
            - `dict[str : list[int]]` Dictionary of dnaA-box indexes on the DNA sequence.
        - `max_mismatches`:
            - `int` Maximum allowed mismatches in the proposed dnaA-boxes.
        - `genes_of_interest`:
            - `list[str]` List of user-given genes that relate to the location of oriCs.
        - `gene_locations`:
            - `list[tuple[str, Peak]]` List of locations of the `genes_of_interest`.
        - `NCBI_oriC`:
            - `list[tuple[str, Peak]]` Only created if the file/accesssion analysed had an annotated 'rep_origin'
        - `windows`:
            - `list[int]` Windows surrounding potential oriCs. Used for filtering, distance approximation, and generalisations
        - `max_point_spread`:
            - `float` Maximum distance between points in a component as a fraction of the length of the sequence
        - `model`:
            - `sklearn.classifier` Scikit-learn classifier used for predicting the significance of the found oriCs.
        - `seq`:
            - `str` String representation of the DNA sequence
        - `seq_len`:
            - `int` Length of the DNA sequence
        - `x`, `y`, `z`, `gc`:
            - `numpy.ndarray` Three components of the Z-curve as well as the GC-skew curve
        - `oriCs`:
            - `list[Peak]` List of Peak object of the predicted origin locations. Consult its scores and prediction for analysis of its validity.
        - `oriC_middles`:
            - `list[int]` List of only the middle posistion point of each predicted oriC in the same order as `oriCs`.
        - `Z_scores`, `G_scores`, `D_scores`:
            - `list[float]` List of scores based on the analysis of the Z-curve, gene locations, and dnaA-box regions. Consult the paper
            for more information. In the same order as `oriCs`.
        - `predictions`:
            - `list[float]` List of decision values decided by the given `model`. In the same order as `oriCs`.

        --------------------------
        Example:
        --------

        >>> # To instantiate the ORCA object call:
        >>> orca = ORCA.from_accession("NC_example.3", "your@email.com")
        >>> # To find oriCs:
        >>> orca.find_oriCs()
        >>> print(orca.accession)
        NC_example.3
        '''

        self.__check_run_properly("accession", "This object has not been instantiated properly.")

        # Step 1: Calculate disparity curves + finding dnaA-box regions
        self.calculate_disparity_curves()

        # Step 2: Finding potential oriCs based on disparity curves
        peaks_of_interest = self.analyse_disparity_curves()
        oriCs, Z_scores = self.calculate_Z_scores(peaks_of_interest)

        # Step 3: DnaA-box analysis
        if len(self.dnaa_dict.keys()) == 0:
            warnings.warn(f'''Accession: {self.accession + "." + str(self.version)}.
                No DnaA-boxes were found: {self.input_dnaa_boxes}
                Allowed for {self.max_mismatches} mismatches.
                This will affect the validity of the prediction!
            ''')
        D_scores = self.calculate_D_scores(oriCs)

        # Step 4: Gene-location analysis
        if len(self.gene_locations) == 0:
            warnings.warn(f'''Accession: {self.accession + "." + str(self.version)}.
                None of the genes of interest were found in the given data: {self.genes_of_interest}
                This will affect the validity of the prediction!
            ''')
        G_scores = self.calculate_G_scores(oriCs)

        # Step 5: Machine Learning model decision function.
        decisions = [None for i in range(len(oriCs))]
        if self.model is not None:
            decisions = self.model.predict_proba(np.asarray(
                [Z_scores, G_scores, D_scores]).T)[:, 1].tolist()
        oriC_middles = [oriC.middle for oriC in oriCs]

        # Step 6: Setting the last variables to the proper values
        for i in range(len(oriCs)):
            oriCs[i].z_score = Z_scores[i]
            oriCs[i].g_score = G_scores[i]
            oriCs[i].d_score = D_scores[i]
            oriCs[i].decision = decisions[i]

        self.oriCs = oriCs
        self.oriC_middles = oriC_middles
        self.Z_scores = Z_scores
        self.G_scores = G_scores
        self.D_scores = D_scores
        self.predictions = decisions
        self._set_best_pot_oriC_idx()

        if show_info:
            self.pretty_print_results()
        if show_plot:
            self.plot_oriC_curves()


def example_use() -> ORCA:
    """Example showcase. See code"""
    import joblib

    email = 'real@email.address'

    # Provided model is compressed due to GitHub's file size limits.
    # Pickle file for the provided model. Model is around 189 MB uncompressed and 30 MB compressed,
    # so it is too large to upload to GitHub.
    model = joblib.load("data/output/machine_learning/ORCA_RFC_model.pkl.gz")

    # orca = ORCA.from_pkl("data/input/NC_000913_3.pkl", model=model)
    orca = ORCA.from_accession("NC_014248", email=email, model=model)
    orca.find_oriCs(show_info=True, show_plot=True)
    return orca


if __name__ == '__main__':
    orca = example_use()
