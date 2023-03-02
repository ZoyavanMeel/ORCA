# Libraries
import os, warnings, pickle, re
from itertools import product, combinations
from typing import Generator

import numpy as np
from sklearn.base import is_classifier
from Bio import Seq, SeqIO

# Self-made modules
from Handlers import CurveHandler, FileHandler
from Peak import Peak
import plotter_functions as pf


class ORCA:
    """
    Private default constructor. Use `from_gbk()`, `from_pkl()`, or `from_accession()` for proper functionality of ORCA.

    Accepted **kwargs:
        - `dnaa_boxes`          : If None, will use the [consensus DnaA-box](https://doi.org/10.1093/bib/bbn031): `TTAT(A|T|C|G)CACA`.
                                Else, provide a list of 9 base strings. See the `get_dnaa_boxes` function in `helper_functions.py` for some more examples of dnaA-boxes.
                                Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
        - `max_mismatches`      : Maximum allowed mismatches allowed in a dnaa_box for it still to be read as such. Recommended: 0; recommended max: 2.
        - `genes_of_interest`   : List of gene names to consider as 'oriC-proximal' and use for helping estimate the location of the oriC.
        - `max_group_spread`    : Maximum spread a group can have when looking for connected groups. Default is 5 % of the total chromosome length
        - `windows`             : The windows around around peaks of skew curves to consider. Defaults are 1, 3 and 5 % of the total chromosome length.
        - `model`               : A fitted scikit-learn classifier. Recommended to use the one provided on [GitHub](https://github.com/ZoyavanMeel/ORCA/).
    """
    def __init__(self, **kwargs):
        user_args = ORCA.setup_variables(**kwargs)
        for key, value in user_args.items():
            setattr(self, key, value)


    @classmethod
    def from_gbk(cls, path: str, **kwargs):
        """
        Instatiates an ORCA object from a GenBank (gbk) file.

        Parameters:
        - `path` : string of path to the GenBank file to analyse.
        - `**kwargs` : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        - `ORCA` : an ORCA object with properly loaded parameters.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} does not exist.")
        orca = cls(**kwargs)
        file_args = FileHandler.parse_SeqRecord(SeqIO.read(path, 'gb'), orca.genes_of_interest)
        for key, value in file_args.items():
            setattr(orca, key, value)
        return orca


    @classmethod
    def from_pkl(cls, path: str, **kwargs):
        """
        Instatiates an ORCA object from a pickled Biopython SeqRecord.
        
        Parameters:
        - `path` : string of path to the pickled SeqRecord of a GenBank file.
        - `**kwargs` : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        - `ORCA` : an ORCA object with properly loaded parameters.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} does not exist.")
        orca = cls(**kwargs)
        with open(path,"rb") as fh:
            record = pickle.load(fh)
            fh.close()
            if not isinstance(record, SeqIO.SeqRecord):
                raise ValueError("Pickled object must be a SeqRecord object")
            file_args = FileHandler.parse_SeqRecord(record, orca.genes_of_interest)
            for key, value in file_args.items():
                setattr(orca, key, value)
            return orca


    @classmethod
    def from_accession(cls, accession: str, email: str, api_key: str = None, **kwargs):
        """
        Instatiates an ORCA object from a given accession number.
        
        Parameters:
        - `accession` : accession to analyse. ORCA will fetch and parse the required file.
        - `email`     : your email address. Always tell NCBI who you are.
        - `api_key`   : Optional. Only necessary if you are doing more than 10 requests per second.
        - `**kwargs`  : Check the ORCA class documentation for valid keyword arguments.

        Returns:
        - `ORCA` : an ORCA object with properly loaded parameters.
        """
        orca = cls(**kwargs)
        if email is None:
            raise ValueError(f"""
                NCBI requires you to specify an email address to make use of their services.
                \tCreate an NCBI account at: https://www.ncbi.nlm.nih.gov/
                \tGenerate an API_key at: https://www.ncbi.nlm.nih.gov/account/settings/ (optional)
            """)
        with FileHandler.fetch_file(accession, email, api_key, rettype='gbwithparts') as handle:
            file_args = FileHandler.parse_SeqRecord(SeqIO.read(handle, 'gb'), orca.genes_of_interest)
            handle.close()
            for key, value in file_args.items():
                setattr(orca, key, value)
            return orca


    @staticmethod
    def setup_variables(**kwargs) -> dict:
        """
        Handle the user input parameters and throw errors and warnings if they are not valid.
        
        Accepted **kwargs:
        - `dnaa_boxes`          : If None, will use the [consensus DnaA-box](https://doi.org/10.1093/bib/bbn031): `TTAT(A|T|C|G)CACA`.
                                Else, provide a list of 9 base strings. See the `get_dnaa_boxes` function in `helper_functions.py` for some more examples of dnaA-boxes.
                                Example input: `['AAAAAAAAA', 'TTTTTTTTT']`.
        - `max_mismatches`      : Maximum allowed mismatches allowed in a dnaa_box for it still to be read as such. Recommended: 0; recommended max: 2.
        - `genes_of_interest`   : List of gene names to consider as 'oriC-proximal' and use for helping estimate the location of the oriC.
        - `max_group_spread`    : Maximum spread a group can have when looking for connected groups. Default is 5 % of the total chromosome length
        - `windows`             : The windows around around peaks of skew curves to consider. Defaults are 1, 3 and 5 % of the total chromosome length.
        - `model`               : A fitted scikit-learn classifier. Recommended to use the one provided on [GitHub](https://github.com/ZoyavanMeel/ORCA/).
        """

        defaults = {
            'dnaa_boxes'        : ['TTATACACA', 'TTATTCACA', 'TTATCCACA', 'TTATGCACA'],
            'genes_of_interest' : ['dnaA', 'dnaN'],
            'windows'           : [0.01, 0.03, 0.05],
            'max_group_spread'  : 0.05,
            'max_mismatches'    : 0,
            'model'             : None,
        }

        # Update defaults with provided values
        args = {key: kwargs[key] if key in kwargs else defaults[key] for key in defaults}

        # Simple type enforcing. Overriding default parameters is okay, but provide empty lists if you do.
        assert args['model'] is None or is_classifier(args['model']), 'model must be a scikit-learn classifier.'
        assert isinstance(args['windows'], list), 'windows must be a list of integers.'
        assert isinstance(args['dnaa_boxes'], list), 'dnaa_boxes must be a list of strings.'
        assert isinstance(args['max_mismatches'], int), 'max_mismatches must be an integer.'
        assert isinstance(args['max_group_spread'], float), 'max_group_spread must be an integer.'
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
        """

        #               Sequences                       DOI                                            Year  Notes
        # consensus_1 = ['TTAT(A|C)CA(A|C)A']         # https://doi.org/10.1016/0092-8674(84)90284-8  (1984) The first consensus
        # consensus_2 = ['TGTG(G|T)ATAAC']            # https://doi.org/10.1016/0022-2836(85)90299-2  (1985) Matsui-box
        # consensus_3 = ['TTAT(A|T|C|G)CACA']         # https://doi.org/10.1093/dnares/dsm017         (2007) In (roughly) all eubacteria, not just B. subtilis
        # consensus_4 = [                             # https://doi.org/10.1093/emboj/16.21.6574      (1997) Affinity study
        #     'TTATCCACA', 'TTTTCCACA', 'TTATCAACA',
        #     'TCATTCACA', 'TTATACACA', 'TTATCCAAA'
        # ]
        # consensus_5 = [                             # https://doi.org/10.1007/BF00273584            (1991) Only in E. coli K12. Do not use.
        #     '(T|C)(T|C)(A|T|C)T(A|C)C(A|G)(A|C|T)(A|C)'
        # ]

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
                raise ValueError(f'\n\tInput string: \'{box}\' contains forbidden characters. Only use the four nucleotide bases: A, T, C, and G.')
            # Also check reverse complement for box on other strand
            boxes.append(str(Seq.Seq(box).reverse_complement()))

        # Get all unique strings while allowing for max. 2 mismatches.
        mismatch_boxes = list(set(boxes))
        for box in boxes:
            for i in range(abs(max_mismatches)):
                mismatch_boxes.extend(list(generate_mismatched_strings(box, i+1)))
        return mismatch_boxes


    def calculate_disparity_curves(self):
        '''
        Z-curve and GC-skew calculation and k-mer indexing. In one function so only one iteration of the sequence is necessary.\n
        Sets:
            `x`, `y`, `z`, `gc` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
            `dnaa_dict`         : Dictionary of middle indexes of dnaa-boxes in `seq`
        '''
        k = len(self.dnaa_boxes[0])
        x, y, z, gc = [], [], [], []
        a, c, t, g  = 0, 0, 0, 0

        raw_dict = {}

        for i in range(self.seq_len):
            base = self.seq[i]
            if base == "A": a +=1
            elif base == "C": c +=1
            elif base == "G": g +=1
            elif base == "T": t +=1

            gc.append(g - c)
            x.append( (a + g) - (c + t) ) # Purine vs Pyrimidine
            y.append( (a + c) - (g + t) ) # Amino vs Keto
            z.append( (a + t) - (c + g) ) # Weak vs Strong Hydrogen Bonds

            kmer = self.seq[i:i+k] if i <= self.seq_len - k else self.seq[i:] + self.seq[:k-(self.seq_len-i)]
            if kmer in self.dnaa_boxes:
                mid = i+k//2+1 if i+k//2+1 <= self.seq_len else i+k//2+1 - self.seq_len
                try: raw_dict[kmer].append(mid)
                except KeyError: raw_dict[kmer] = [mid]

        (self.x, self.y, self.z, self.gc, self.dnaa_dict) = (np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(gc), raw_dict)


    def analyse_disparity_curves(self) -> list[Peak]:
        '''
        Analyse the three disparity curves related to finding the oriC of the sequence.
        The Z-curve can only be analysed after calling `calc_disparities()`.
        '''
        peaks = []
        for fraction in self.windows:
            window_size = int(self.seq_len * fraction)
            peaks_x  = CurveHandler.process_curve(self.x, 'min', window_size=window_size)
            peaks_y  = CurveHandler.process_curve(self.y, 'max', window_size=window_size)
            peaks_gc = CurveHandler.process_curve(self.gc, 'min', window_size=window_size)
            peaks.extend( [j for i in CurveHandler.curve_combinations( (peaks_x, peaks_y, peaks_gc), self.seq_len ) for j in i] )
        return peaks


    def calculate_Z_scores(self, peaks: list[Peak]) -> tuple[list[Peak], list[int]]:
        '''Finding connected components in undirected graph with a depth-first search to merge Z-curve oriCs'''

        adj_mat = Peak.get_adjacency_matrix(peaks)
        connected_groups = Peak.get_connected_groups(peaks, adj_mat, int(self.seq_len*self.max_group_spread))
        connected_groups.sort(key=lambda x:len(x), reverse=True)

        total_pot_oriCs = len( [y for x in connected_groups for y in x] )
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
        # LOOK INTO THIS: what if all genes of interest are far apart?
        if len(self.gene_locations) == 0:
            return [0] * len(peaks)
        gene_peaks = [name_loc_tuple[1] for name_loc_tuple in self.gene_locations]
        matrix = Peak.get_adjacency_matrix(peaks, gene_peaks)
        if np.min(matrix) == np.max(matrix):
            return [1] * matrix.shape[1]
        norm_mat = (matrix - np.min(matrix)) / (np.max(matrix) - np.min(matrix))
        return [1 - x for x in np.mean(norm_mat, axis=1)]


    def find_oriCs(self, show_info: bool = False, show_plot: bool = False):
        '''
        Locates potential oriCs on circular bacterial chromosomes based on Z-curve and GC-skew analysis, dnaA box analysis, and dnaA/dnaN gene locations.
        Three default window_sizes are used: 1, 3 and 5 % of the total genome length. See the README-file in the [GitHub repository](https://github.com/ZoyavanMeel/ORCA/)
        for more information or consult ORCA.pdf for more extensive results and analyses of ORCA.

        Parameters:
        - `show_info`           : If True, prints info of ALL found oriCs. Good and bad.
        - `show_plot`           : If True, shows plot of ALL found oriCs. Good and bad. Should not be used for analysis -> Make a separate plot for the best oriCs.

        Return:
        - `properties`   : Dictionary with properties of all oriC-like regions.
                           NOTE: oriCs are NOT sorted by importance. Recommended way to rank: learning machine decision.
        - `'oriCs`       : List of Peak-objects. Each Peak has a index position on the given sequence and scores based on the analyses (see: ORCA.pdf).
        - `'dnaA_boxes'` : Dictionary with dnaA-box 9-mers as keys and lists of position indices on the given DNA as values.
                           The indices refer to the position of the 5th base in the 9-mer.
        - etc.
        '''

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
                Will not use DnaA-boxes in prediction.
            ''')
        D_scores = self.calculate_D_scores(oriCs)

        # Step 4: Gene-location analysis
        if len(self.gene_locations) == 0:
            warnings.warn(f'''Accession: {self.accession + "." + str(self.version)}.
                None of the genes of interest were found in the given data: {self.genes_of_interest}
                Will not use gene locations in prediction.
            ''')
        G_scores = self.calculate_G_scores(oriCs)

        # Step 4: Machine Learning model decision function.
        decisions = [None for i in range(len(oriCs))]
        if self.model is not None:
            decisions = self.model.decision_function(np.asarray([Z_scores, G_scores, D_scores]).T).tolist()
        oriC_middles = [oriC.middle for oriC in oriCs]

        # Step 5: Setting the last variables to the proper values
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

        if show_info:
            print('Predictions :', decisions)
            print('Z-scores    :', Z_scores)
            print('G-scores    :', G_scores)
            print('D-scores    :', D_scores)
            print('oriCs       :', oriC_middles)
        if show_plot:
            pf.plot_Z_curve_2D([self.x, self.y, self.gc], oriC_middles, ['$x_n$', '$y_n$', '$g_n$'])


if __name__ == '__main__':
    import joblib

    acc = 'NC_000913'
    email = 'no_need_for_a_real@email_address.com'
    model = joblib.load("Machine_learning/75_train_model.pkl")

    orca = ORCA.from_pkl("Test/NC_000913_3.pkl", model=model)
    orca.find_oriCs(True, True)

    # orca_dict = ORCA.find_oriCs(
    #     accession='NC_000913', # E. coli K-12       #'NC_000117'
    #     email=email,
    #     api_key=None,
    #     model=model,
    #     show_plot=True,
    #     show_info=True
    # )

    # for key in orca_dict.keys():
    #     if key != 'z_curve' and key != 'dnaA_boxes':
    #         print(key + ':\t', orca_dict[key])