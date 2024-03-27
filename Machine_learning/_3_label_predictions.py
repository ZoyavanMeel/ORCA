"""This script is made in order to label the values of the ORCA scores according to the ground truths (DoriC or experimentally verified)."""

from Peak import Peak
from _1_download_dataset import load_data
import os
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


GROUND_TRUTHS = ['data/input/DoriC_chromosome_circular.csv', 'data/input/DoriC_complete_circular.csv']
INDICATOR_VALUES_CSV = "data/output/doric_set_no_model_orca.csv"
OUT_FILE_LABELS = "data/output/machine_learning/labels.csv"


def label(indicator_values_csv: str, ground_truth_csvs: list[str], output_csv: str) -> None:
    """Compare all possible found oriCs by ORCA to the ground truth oriCs for the corresponding accession.
    If the predicted origin is further than 2.5% of the total genome length removed from the ground truth,
    it gets labelled as a `False` prediction. If it is within this distance, it gets labelled as a `True` prediction.

    This method outputs a CSV with the given name that can be used for training ML models."""

    gt_df = load_data(*ground_truth_csvs)
    iv_df = pd.read_csv(indicator_values_csv)

    gt_df.drop(labels=['oric_seq', 'lineage'], axis=1, inplace=True)
    iv_df['nc'] = iv_df['accession'] + '.' + iv_df['version'].map(str)

    iv_df['nc'].astype('string')
    gt_df['nc'].astype('string')

    gt_df = gt_df.merge(iv_df, on='nc', how='right')[gt_df.columns.to_list() + ['seq_len']]

    gt_df['middles'] = gt_df.apply(lambda x: Peak.calc_middle(x.oric_start, x.oric_end, x.seq_len), axis=1)

    label_dict = {
        'accession': [],
        'version': [],
        # e.g. if `total_pot` = 5, then pot_oriC_num is unique (1 through 5) for each of those 5 potential oriCs
        'pot_oriC_num': [],
        'Z_score': [],
        'G_score': [],
        'D_score': [],
        'total_pot': [],  # how many other potential oriCs did ORCA think were on this genome?
        'correct': []
    }

    for _, sample_original in iv_df.iterrows():
        sample = sample_original.copy()
        sample.dropna(inplace=True)
        ORCA_oriC_cols = [i for i in sample.axes[0] if 'oriC_middle_' in i]
        ORCA_peaks = [Peak(int(sample[ORCA_oriC_cols[i]]), sample['seq_len'], 0) for i in range(len(ORCA_oriC_cols))]

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
            gt_middles = gt_df.loc[gt_df['nc'] == sample['nc']]['middles']

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
    label(INDICATOR_VALUES_CSV, GROUND_TRUTHS, OUT_FILE_LABELS)


if __name__ == '__main__':
    main()
