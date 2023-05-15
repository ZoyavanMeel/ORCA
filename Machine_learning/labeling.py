"""This script is made in order to label the values of the ORCA scores according to the ground truths."""
import os, sys
import pandas as pd

from predict_dataset_locally import predict_dataset

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from Peak import Peak

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

INPUT_PATH_ACC       = 'data\input\dataset.csv'   # Only used for getting the accession numbers.
INPUT_PATH_GT        = 'data\input\dataset.csv'   # Ground-Truth values
CSV_OUT_FOLDER       = 'data\output\experimental_set_predictions_orca'
OUT_FILE_PREDICTIONS = "data\output\experimental_set_predictions_orca\orca.csv"
OUT_FILE_LABELS      = "data\output\machine_learning\labels.csv"
PARALLEL             = True                 # Whether to do parallel processing for the prediction step.
CPUS                 = os.cpu_count() - 1   # Amount of CPUs to use if using parallel processing.
MAX_ORICS            = 10                   # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'no_need_for_a_real@email_address.com'
API_KEY = None
MODEL   = None


def label(indicator_values_csv: str, ground_truth_csv: str, output_csv: str) -> None:
    """Compare all possible found oriCs by ORCA to the ground truth oriCs for the corresponding accession.
    If the predicted origin is further than 2.5% of the total genome length removed from the ground truth,
    it gets labelled as a `False` prediction. If it is within this distance, it gets labelled as a `True` prediction.
    
    This method outputs a CSV with the given name that can be used for training ML models."""

    gt_df = pd.read_csv(ground_truth_csv)
    iv_df = pd.read_csv(indicator_values_csv)

    gt_df['middles'] = gt_df.apply(lambda x: Peak.calc_middle(x.begin, x.end, x.seq_len), axis=1)

    label_dict = {
        'accession'    : [],
        'version'      : [],
        'pot_oriC_num' : [],
        'Z_score'      : [],
        'G_score'      : [],
        'D_score'      : [],
        'total_pot'    : [], # how many other potential oriCs did ORCA think were on this genome?
        'correct'      : []
    }

    for _, sample_original in iv_df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)
        ORCA_oriC_cols = [i for i in sample.axes[0] if 'oriC_middle_' in i]
        ORCA_peaks = [ Peak( int(sample[ORCA_oriC_cols[i]]), sample['seq_len'], 0) for i in range(len(ORCA_oriC_cols)) ]

        # Loop over each potential oriC that ORCA found to check if it is correct
        for i, peak in enumerate(ORCA_peaks):

            label_dict['accession'].append(sample.accession)
            label_dict['version'].append(sample.version)
            label_dict['pot_oriC_num'].append(i+1)
            label_dict['Z_score'].append(sample[f'Z_score_{i}'])
            label_dict['G_score'].append(sample[f'G_score_{i}'])
            label_dict['D_score'].append(sample[f'D_score_{i}'])
            label_dict['total_pot'].append(len(ORCA_peaks))

            # Ground truth middle value. One accession can have multiple oriCs (e.g. bipartite oriC)
            gt_middles = gt_df.loc[gt_df['accession'] == sample.accession]['middles']

            # Check correctness of each potential oriC and label them as such
            # Loop over each ground truth oriC in the accession (usually 1) to check if the potential oriC corresponds to any of them
            correctness = False
            for _, gt_middle in gt_middles.items():
                if peak.calc_dist(gt_middle) < sample.seq_len * 0.025:
                    correctness = True
            label_dict['correct'].append(correctness)

    # Output the CSV
    # print(pd.DataFrame(label_dict))
    pd.DataFrame(label_dict).to_csv(output_csv, index=False)
    


def main():
    predict_dataset(INPUT_PATH_ACC, CSV_OUT_FOLDER, OUT_FILE_PREDICTIONS, PARALLEL, CPUS, MAX_ORICS, EMAIL, API_KEY, MODEL)
    label(OUT_FILE_PREDICTIONS, INPUT_PATH_GT, OUT_FILE_LABELS)


if __name__ == '__main__':
    main()