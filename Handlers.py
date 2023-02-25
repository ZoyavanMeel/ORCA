import os, shutil, csv, re
from itertools import combinations, product
from typing import TextIO, Union, Tuple, Generator
from urllib.error import HTTPError, URLError

import numpy as np
from Bio import Seq, SeqIO, Entrez

from Curve import Curve
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
            (api_message_1, api_message_2) = (f' with API_key \'{api_key}\'', ' and if your API_key if correctly typed and linked to the given email') if api_key is not None else ('', '')
            raise ValueError(f'Unable to fetch accession: \'{accession}\' using \'{email}\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.') from BadRequest
        except URLError as NoConnection:
            raise ConnectionError('You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.') from NoConnection


    @staticmethod
    def read_FASTA(handle: Union[TextIO, str]) -> Tuple[str, str]:
        """Read a FASTA file and returns the name and sequence of only the first entry in the file"""
        Seq_records = SeqIO.parse(handle, 'fasta')
        Seq_obj = next(Seq_records)
        return Seq_obj.id, Seq_obj.seq


    @staticmethod
    def read_gene_info(handle: TextIO, genes_list: list[str]) -> dict:
        """
        Read FASTA-file acquired with rettype='fasta_cds_na'. `genes_list` is a list of names of genes to extract.\n
        Return:
        - The features of the genes in the `genes_list`.
        - Total number of genes in the `handle`.
        """

        obj = SeqIO.parse(handle, 'fasta')
        genes = [gene.lower() for gene in genes_list]
        genes_dict = {}
        for gene in obj:
            features = [x.split('=') for x in re.findall(r"\[(.*?)\]", gene.description) if '=' in x]
            feature_dict = {feature[0] : feature[1] for feature in features}
            try: gene_name = feature_dict.pop('gene')
            except KeyError: continue # be careful about trying to get an alternative name
            if gene_name.lower() in genes:
                genes_dict.update({gene_name : feature_dict})
        return genes_dict


    @staticmethod
    def move_fastas(db_loc, on_cluster=True, split=4):
        '''
        Split one folder into 'split' amount of folders with roughly the same amount of files in them.
        Can be used for easier parallel processing. Instead of one big job with 4000 samples. Run 4 jobs with 1000 samples at the same time.
        '''
        if on_cluster:
            path = db_loc + '/bacteria'
            samples = os.listdir(path)
        else:
            path = db_loc + '/chromosomes_only'
            samples = os.listdir(path)

        samples_per_split = len(samples)//split
        for i in range(split):
            new_folder = db_loc + '/' + str(i)
            os.mkdir(new_folder)
            start = 0 + i*samples_per_split
            stop  = (i+1)*samples_per_split if i != split-1 else None
            for sample in samples[start:stop]:
                shutil.move(path + '/' + sample, new_folder + '/' + sample)


    @staticmethod
    def download(accession: str, output_folder: str, email: str = None, api_key: str = None):
        '''Download the proper file types for large dataset analysis into the output_folder.'''
        Entrez.email = email
        if api_key is not None:
            Entrez.api_key = api_key

        if not os.path.exists(output_folder + f'/gene_info_files/{accession}.fasta'):
            check_handle = Entrez.efetch(db='nuccore', id=accession, rettype='gb', retmode='text')
            annotations = next(SeqIO.parse(check_handle, 'genbank').records).annotations
            if annotations['taxonomy'][0].lower() == 'bacteria' and annotations['topology'].lower() == 'circular':
                handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta_cds_na", retmode="text")
                with open(output_folder + f'/gene_info_files/{accession}.fasta', 'w') as fh:
                    lines = handle.read()
                    fh.write(lines)
                    handle.close()
                    fh.close()
                handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
                with open(output_folder + f'/bacteria/{accession}.fasta', 'w') as fh:
                    lines = handle.read()
                    fh.write(lines)
                    handle.close()
                    fh.close()


    @staticmethod
    def merge_csvs(file_folder: str, merged_csv: str, fieldnames: list[str], length: int = -1, headers: bool = False):
        '''
        Merge multiple csvs into one csv.
        Arguments:
            file_folder : path to folder with csvs that have to be merged
            merged_csv  : name of the merged csv
            fieldnames  : list of fieldnames for the merged_csv
            length      : amount of rows in each single csv. -1, if the length of each single
                        csv is not known or not the same.
            headers     : if the single csvs have headers or not
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



class SequenceHandler:
    '''
    Sequence attributes:
    - `length` : length of the sequence
    - `gc_conc` : (G-count + C-count) / length
    - `pot__boxes` : list of all *possible* DnaA binding sites
    - `x`, `y`, `z` : all components of the Z-curve obtained from the sequence as `Curve` objects
    - `gc` : GC-skew of the sequence
    - `dnaa_dict` : dictionary of all DnaA-boxes actually in the sequence
    
    The DNA sequence string itself is not saved.
    '''
    def __init__(self, args: dict):
        # Sequence fetching and reading
        if args['genome_fasta'] is None:
            seq_handle = FileHandler.fetch_file(args['accession'], args['email'], args['api_key'], 'fasta')
        else:
            seq_handle = args['genome_fasta']
        args['accession'], sequence = FileHandler.read_FASTA(seq_handle)

        self.length  = len(sequence)
        self.windows = args['windows']
        self.gc_conc = ( sequence.count('G') + sequence.count('C') ) / self.length
        self.max_group_spread = args['max_group_spread']

        # Parsing sequence properties
        self.pot_boxes = self.get_dnaa_boxes(box_list=args['dnaa_boxes'], max_mismatches=args['max_mismatches'])
        self.x, self.y, self.z, self.gc, self.dnaa_dict = self.calc_disparities(sequence)


    def calc_disparities(self, seq: str) -> tuple[Curve, Curve, Curve, Curve, dict]:
        '''
        Z-curve and GC-skew calculation and k-mer indexing. In one function so only one iteration of the sequence is necessary.\n
        Parameters:
            `seq`        : string DNA sequence
            `k`          : length of k-mer. Should be same length as every dnaa-box
            `dnaa_boxes` : set of dnaa-box regions.\n
        Return:
            `x`, `y`, `z`, `gc` : 1D-np.arrays of the three Z-curve components and 1D-np.array of the GC-skew
            `dnaa_dict`         : Dictionary of starting indexes of dnaa-boxes in `seq`
        '''
        k = len(list(self.pot_boxes)[0])
        x, y, z, gc = [], [], [], []
        a, c, t, g  = 0, 0, 0, 0

        raw_dict = {str: list}
        seq_len  = len(seq)

        for i in range(seq_len):
            base = seq[i]
            if base == "A": a +=1
            elif base == "C": c +=1
            elif base == "G": g +=1
            elif base == "T": t +=1

            gc.append(g - c)
            x.append( (a + g) - (c + t) ) # Purine vs Pyrimidine
            y.append( (a + c) - (g + t) ) # Amino vs Keto
            z.append( (a + t) - (c + g) ) # Weak vs Strong Hydrogen Bonds

            kmer = seq[i:i+k] if i <= seq_len - k else seq[i:] + seq[:k-(seq_len-i)]
            mid = i+k//2+1 if i+k//2+1 <= seq_len else i+k//2+1 - seq_len
            try: raw_dict[kmer].append(mid)
            except KeyError: raw_dict[kmer] = [mid]

        # & assumes .keys() as a set (which it should as dict keys are unique). .intersection() assumes .keys as a list and set.intersection(list) has a worse time-complexity. https://wiki.python.org/moin/TimeComplexity
        keys = self.pot_boxes & raw_dict.keys()
        dnaa_dict = {key : raw_dict[key] for key in keys}
        del raw_dict
        return Curve(np.asarray(x), 'min'), Curve(np.asarray(y), 'max'), Curve(np.asarray(z), ''), Curve(np.asarray(gc), 'min'), dnaa_dict


    def analyse_Z_curve(self) -> list[Peak]:
        '''Analyse the three disparity curves related to finding the oriC of the sequence.'''        
        peaks = []
        for fraction in self.windows:
            window_size = int(self.length * fraction)
            peaks_x  = self.x.process_array(window_size=window_size)
            peaks_y  = self.y.process_array(window_size=window_size)
            peaks_gc = self.gc.process_array(window_size=window_size)
            peaks.extend( [j for i in self.curve_combinations( (peaks_x, peaks_y, peaks_gc) ) for j in i] )
        return peaks


    def __merge_oriCs(self, groups: list) -> tuple[list[Peak], list[int]]:
        '''Finds the average index of a group and returns those values. groups is a nested-list'''
        mutable = sorted( groups, key=lambda x:len(x), reverse=True )
        total_pot_oriCs = len( [y for x in mutable for y in x] )
        oriCs, Z_scores = [], []
        window_size = int(self.length*self.windows[-1])

        group: list[Peak]
        for group in mutable:
            group.sort()
            for i in range(len(group)):
                if (group[-1] - group[i]).middle >= (self.length-1)/2:
                    group[i].middle += self.length-1
            avg_val = sum(group)//len(group)
            if avg_val > self.length-1:
                avg_val -= self.length-1
            oriCs.append(Peak(avg_val, self.length, window_size))
            Z_scores.append(len(group)/total_pot_oriCs)
        return oriCs, Z_scores


    def calculate_Z_scores(self, peaks: list[Peak]) -> tuple[list[Peak], list[int]]:
        ## Finding connected components in undirected graph with a depth-first search to merge Z-curve oriCs
        matrix_pot_oriCs = Peak.get_adjacency_matrix(peaks)
        connected_groups = Peak.get_connected_groups(peaks, matrix_pot_oriCs, int(self.length*self.max_group_spread))
        return self.__merge_oriCs(connected_groups)


    def calculate_D_scores(self, peaks: list[Peak]) -> list[Peak]:
        '''Process the location of dnaa_boxes and rank potential oriCs based on most surrounding dnaa_boxes.'''
        if len(self.dnaa_dict.keys()) != 0:
            contains_boxes = []
            all_pos = [pos for pos_list in self.dnaa_dict.values() for pos in pos_list]
            for oriC in peaks:
                count = 0
                for pos in all_pos:
                    if oriC.contains_point(pos):
                        count += 1
                contains_boxes.append(count)
            return [x/sum(contains_boxes) if sum(contains_boxes) != 0 else 0 for x in contains_boxes], False
        else:
            return [0] * len(peaks), True


    def curve_combinations(self, peaks_list: tuple[list["Peak"]]) -> list:
        '''Get every matched_peaks combination for x, y, and gc.'''
        oriC_locations_list = []
        for peaks_i, peaks_j in combinations(peaks_list, 2):
            matched_peaks  = Peak.match_peaks(peaks_i, peaks_j)
            oriC_locations_list.append(
                [Peak( Peak.get_middle(matches[0], matches[1]), self.length, peaks_list[0][0].window_size ) for matches in matched_peaks]
            )
        return oriC_locations_list


    @staticmethod
    def get_dnaa_boxes(box_list: list, max_mismatches: int = 2) -> set:
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

        #             Sequences                       DOI                                            Year  Notes
        consensus_1 = ['TTAT(A|C)CA(A|C)A']         # https://doi.org/10.1016/0092-8674(84)90284-8  (1984) The first consensus
        consensus_2 = ['TGTG(G|T)ATAAC']            # https://doi.org/10.1016/0022-2836(85)90299-2  (1985) Matsui-box
        consensus_3 = ['TTAT(A|T|C|G)CACA']         # https://doi.org/10.1093/dnares/dsm017         (2007) In (roughly) all eubacteria, not just B. subtilis
        consensus_4 = [                             # https://doi.org/10.1093/emboj/16.21.6574      (1997) Affinity study
            'TTATCCACA', 'TTTTCCACA', 'TTATCAACA',
            'TCATTCACA', 'TTATACACA', 'TTATCCAAA'
        ]
        consensus_5 = [                             # https://doi.org/10.1007/BF00273584            (1991) Only in E. coli K12. Do not use.
            '(T|C)(T|C)(A|T|C)T(A|C)C(A|G)(A|C|T)(A|C)'
        ]

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
            if len(box) != 9:
                raise ValueError(f'Input string \'{box}\' not of length 9.')
            # Also check reverse complement for box on other strand
            boxes.append(str(Seq.Seq(box).reverse_complement()))

        # Get all unique strings while allowing for max. 2 mismatches.
        mismatch_boxes = list(set(boxes))
        for box in boxes:
            for i in range(abs(max_mismatches)):
                mismatch_boxes.extend(list(generate_mismatched_strings(box, i+1)))
        return set(mismatch_boxes)



class GeneHandler:
    def __init__(self, args: dict):

        # Gene info fetching and reading
        if args['genes_fasta'] is None:
            gene_handle = FileHandler.fetch_file(args['accession'], args['email'], args['api_key'], 'fasta_cds_na')
        else:
            gene_handle = args['genes_fasta']
        self.genes_dict = FileHandler.read_gene_info(gene_handle, args['genes_of_interest'])
        self.length = args['length']
    

    def analyse_gene_locations(self, oriCs: list[Peak]) -> list[Peak]:
        """Returns list of Peaks of the middle position of every gene in the dictionary. `genes_dict` is assumed to be in the format provided by `read_gene_info()`."""
        locations = []
        if len(self.genes_dict.keys()) == 0:
            return locations
        for gene_dict in self.genes_dict.values():
            clean_locs = GeneHandler.handle_location(gene_dict['location'], self.length)
            locations.extend(clean_locs)
        middles = [Peak(Peak.get_middle(loc[0], loc[1], self.length), self.length, 0) for loc in locations]
        return middles
    

    @staticmethod
    def calculate_G_scores(peaks: list[Peak], gene_locations: list[Peak]) -> tuple[list[int], bool]:
        '''Process the location of the genes of interest and rank the potential oriCs based on how close they are to these genes'''
        # LOOK INTO THIS: what if all genes of interest are far apart?
        if len(gene_locations) == 0:
            return [0] * len(peaks), True
        matrix = Peak.get_adjacency_matrix(peaks, gene_locations)
        if np.min(matrix) == np.max(matrix):
            return [1] * matrix.shape[1]
        norm_mat = (matrix - np.min(matrix)) / (np.max(matrix) - np.min(matrix))
        return [1 - x for x in np.mean(norm_mat, axis=1)], False


    @staticmethod
    def handle_location(location: str, seq_len: int) -> list:
        """
        Gene locations come in a billion flavours due to different annotation options. This function handles the vast majority of cases.
        Biopython can parse all location formats (AFAIK), but I can't find how to access their parsers, since they're hidden/private in some of the classes
        """
        handled = []
        if re.search(r'[^joincmplemt<>,\s\d\(\)\.]', location) is not None:
            # LOOK INTO THIS.
            raise ValueError("Location format not supported: " + location)
        loc_groups = GeneHandler.split_location(location)
        for loc_lst in loc_groups:
            if len(loc_lst) == 1:
                loc_coords = re.findall(r'\d+', loc_lst[0])
                if len(loc_coords) == 0:
                    return None
                if len(loc_coords) == 1:
                    loc_coords.append(loc_coords[0])
                loc_coords = [int(x) for x in loc_coords]
                handled.append(loc_coords)
            else: # len(i) > 1
                all_nums = [int(y) for x in loc_lst for y in re.findall(r'\d+', x)]
                relevant_idxes = np.unravel_index(Peak.get_adjacency_matrix(all_nums, seq_len=seq_len).argmax(), (len(all_nums), len(all_nums)))
                first, second = min(all_nums[relevant_idxes[0]], all_nums[relevant_idxes[1]]), max(all_nums[relevant_idxes[0]], all_nums[relevant_idxes[1]])
                handled.append( [first, second] )
        return handled


    @staticmethod
    def split_location(location):
        """
        Splits a `location` string into its location groups.

        (Extreme) e.g.:
        - `location = 'join(complement(1..>2),3..4,<5..6,500),100..200,join(500..600,complement(0..4))'`\n
        returns:
        - `[['1..>2', '3..4', '5..6', '500'], ['100..200'], ['500..600', '0..4']]`
        """
        raw_locs = re.split(r',\s?', location)
        locs = []
        done, i = False, 0
        while not done: # join(...) groups suck
            if i > len(raw_locs) - 1:
                break
            if raw_locs[i].count('(') != raw_locs[i].count(')'):
                group_done, j = False, i+1
                while not group_done:
                    if j > len(raw_locs) - 1:
                        break
                    (group_done, j) = (True, j) if raw_locs[j].count('(') != raw_locs[j].count(')') else (False, j+1)
                locs.append( ','.join([raw_locs[x] for x in range(i, j+1)]) )
                i = j+1
            else:
                locs.append(raw_locs[i])
                i += 1
        locs = [re.findall(r'\d+(?:\.\.[<>]?\d+)?', entry) for entry in locs]
        return locs
