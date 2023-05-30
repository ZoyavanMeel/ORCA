"""This script is made in order to label the oriCs of DoriC according to the ground truths."""
import os, sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
from Peak import Peak

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def label(doric_values_csv: str, ground_truth_csv: str, output_csv: str) -> None:
    """Compare all possible found oriCs by ORCA to the ground truth oriCs for the corresponding accession.
    If the predicted origin is further than 2.5% of the total genome length removed from the ground truth,
    it gets labelled as a `False` prediction. If it is within this distance, it gets labelled as a `True` prediction.
    
    This method outputs a CSV with the given name that can be used for training ML models."""

    gt_df = pd.read_csv(ground_truth_csv)
    dv_df = pd.read_csv(doric_values_csv)

    gt_df['middles'] = gt_df.apply(lambda x: Peak.calc_middle(x.begin, x.end, x.seq_len), axis=1)

    label_dict = {
        'accession'    : [],
        'version'      : [],
        'correct'      : []
    }
    count = 0
    for _, sample in dv_df.iterrows():
        accession, version = tuple(sample.nc.split("."))
        # print(accession, version)
        gt_middles = gt_df.loc[gt_df['accession'] == accession]['middles']
        temp = list(gt_df.loc[gt_df['accession'] == accession]['seq_len'])
        if len(temp) != 0:
            label_dict['accession'].append(accession)
            label_dict['version'].append(version)
            gt_seq_len = temp[0]
            dv_middle = Peak.calc_middle(sample.oric_start, sample.oric_end, gt_seq_len)

            correctness = False
            for _, gt_middle in gt_middles.items():
                if Peak.calc_dist_points(dv_middle, gt_middle, gt_seq_len) < gt_seq_len * 0.025:
                    correctness = True
            label_dict['correct'].append(correctness)
        else:
            count += 1

    # Output the CSV
    print(count)
    print(pd.DataFrame(label_dict))
    # pd.DataFrame(label_dict).to_csv(output_csv, index=False)


if __name__ == '__main__':
    label("data\input\DoriC.csv", "data\input\dataset.csv", "")