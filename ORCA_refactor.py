# Libraries
import warnings, joblib
import numpy as np

# Self-made modules
from Handlers import SequenceHandler, GeneHandler
import plotter_functions as pf


def check_variables(**kwargs) -> dict:
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


def find_oriCs(**kwargs):
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
    args = check_variables(**kwargs)

    # Step 1: Finding potential oriCs based on Z-curve
    seq_handler = SequenceHandler(args)
    peaks_of_interest = seq_handler.analyse_Z_curve()
    oriCs, Z_scores = seq_handler.calculate_Z_scores(peaks_of_interest)

    # Step 2: DnaA-box analysis
    D_scores, warning = seq_handler.calculate_D_scores(oriCs)
    if warning:
        warnings.warn(
            f'Accession: {args["accession"]}.\n' +
            f'No DnaA-boxes were found: {seq_handler.pot_boxes}' +
            'Will not use DnaA-boxes in prediction.\n'
        )

    # Step 3: Gene info fetching and reading
    args.update({'length': seq_handler.length})
    gene_handler = GeneHandler(args)
    try:
        # Gene-location analysis
        gene_locations = gene_handler.analyse_gene_locations(oriCs)
        G_scores, warning = gene_handler.calculate_G_scores(oriCs, gene_locations)
        if warning:
            warnings.warn(
                f'Accession: {args["accession"]}.\n' +
                f'\nNone of the genes of interest were found in the \'genes_fasta\': {args["genes_of_interest"]}' +
                '\nWill not use gene locations in prediction.\n'
            )
    except ValueError as e:
        # Can throw an error in the gene location parsing
        print(e)
        return None
    
    # Step 4: Machine Learning model decision function.
    decisions = [None for i in range(len(oriCs))]
    if args['model'] is not None:
        decisions = args['model'].decision_function(np.asarray([Z_scores, G_scores, D_scores]).T).tolist()
    oriC_middles = [oriC.middle for oriC in oriCs]

    # Setting Z-, G-, and D-scores for each Peak object
    for i in range(len(oriCs)):
        oriCs[i].z_score = Z_scores[i]
        oriCs[i].g_score = G_scores[i]
        oriCs[i].d_score = D_scores[i]
        oriCs[i].decision = decisions[i]

    args.update({
        'oriCs'        : oriCs,
        'oriC_middles' : oriC_middles,
        'Z_scores'     : Z_scores,
        'G_scores'     : G_scores,
        'D_scores'     : D_scores,
        'predictions'  : decisions,
        'z_curve'      : (seq_handler.x.curve, seq_handler.y.curve, seq_handler.z.curve),
        'gc_skew'      : seq_handler.gc.curve,
        'gc_conc'      : seq_handler.gc_conc,
        'dnaA_boxes'   : seq_handler.dnaa_dict
    })

    if args['show_info']:
        print('Predictions :', decisions)
        print('Z-scores    :', Z_scores)
        print('G-scores    :', G_scores)
        print('D-scores    :', D_scores)
        print('oriCs       :', oriC_middles)
    if args['show_plot']:
        pf.plot_Z_curve_2D([seq_handler.x.curve, seq_handler.y.curve, seq_handler.gc.curve], oriC_middles, ['$x_n$', '$y_n$', '$g_n$'])
    return args


if __name__ == '__main__':
    email = 'no_need_for_a_real@email_address.com'
    model = joblib.load('Machine_Learning/75_train_model.pkl')

    orca_dict = find_oriCs(
        accession='NC_000913', # E. coli K-12       #'NC_000117'
        email=email,
        api_key=None,
        model=model,
        show_plot=True,
        show_info=True
    )

    for key in orca_dict.keys():
        if key != 'z_curve' and key != 'dnaA_boxes':
            print(key + ':\t', orca_dict[key])