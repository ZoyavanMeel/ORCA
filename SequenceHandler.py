import re
import numpy as np

from itertools import combinations, product
from Bio import Seq
from typing import Generator

from FileHandler import FileHandler
from Curve import Curve
from Peak import Peak

class SequenceHandler:
    '''
    Sequence attributes:
    - `length` : length of the sequence
    - `gc_conc` : (G-count + C-count) / length
    - `pot__boxes` : list of all *possible* DnaA binding sites
    - `x`, `y`, `z` : all components of the Z-curve obtained from the sequence as 1-D numpy arrays
    - `gc` : GC-skew of the sequence
    - `dnaa_dict` : dictionary of all DnaA-boxes actually in the sequence
    
    The DNA sequence string itself is not saved.
    '''
    
    def __init__(self, args: dict):
        # Sequence fetching and reading
        seq_handle = FileHandler.fetch_file(args['accession'], args['email'], args['api_key'], 'fasta') if args['genome_fasta'] is None else args['genome_fasta']
        args['accession'], sequence = FileHandler.read_FASTA(seq_handle)

        self.length = len(sequence)
        self.gc_conc = ( sequence.count('G') + sequence.count('C') ) / self.length

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
        return Curve(np.asarray(x), 'min', ), Curve(np.asarray(y), 'max'), Curve(np.asarray(z), ''), Curve(np.asarray(gc), 'min'), dnaa_dict


    def curve_combinations(self, peaks_list: tuple[list["Peak"]]) -> list:
        '''Get every matched_peaks combination for x, y, and gc.'''
        oriC_locations_list = []
        for peaks_i, peaks_j in combinations(peaks_list, 2):
            matched_peaks  = Peak.match_peaks(peaks_i, peaks_j)
            oriC_locations_list.append( [Peak(Peak.get_middle(matches[0], matches[1]), self.length, peaks_list[0][0].window_size) for matches in matched_peaks] )
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
                mismatch_boxes.extend(list(SequenceHandler._generate_mismatched_strings(box, i+1)))
        return set(mismatch_boxes)


    @staticmethod
    def _generate_mismatched_strings(string: str, mismatches: int = 2) -> Generator:
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


