"""
Now that we have a means of training a model, we can check ORCA's 'actual' precision and recall

We test it two ways:
- First way: (DoriC_full)
    1. Train a model on all of the DoriC dataset except the experimental data.
    2. Predict the experimental set and only keep the potential oriCs that the model claimed to be true origins.
    3. Check precision and recall of those numbers.

- Second way: (DoriC_70)
    1. Train a model on 70 % of the DoriC dataset without the experimentally determined oigins.
    2. Predict the 30 % test set and the experimental set and only keep the potential oriCs that the model claimed to be true origins.
    3. Check precision and recall of those numbers.
"""

import BioFile as bf
from Peak import Peak
from _4_validate_train_model import load_data_labels_from_df, load_data_labels_from_path
from _1_download_dataset import load_data
import os
import sys
import warnings
import multiprocessing as mp
from functools import partial
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

warnings.filterwarnings(action='ignore', category=UserWarning)

DATA_PATH_DORIC = "data/output/machine_learning/labels.csv"
DATA_PATH_EXP_SET = "data/output/machine_learning/labels_against_NCBI.csv"
MODEL_PATH = 'data/output/machine_learning/DoriC_E_model.pkl'
RANDOM_STATE = 42


def run_DoriC_full() -> None:
    if not os.path.exists(MODEL_PATH):
        DoriC_full(DATA_PATH_DORIC, DATA_PATH_EXP_SET)
    RFC_model = joblib.load(MODEL_PATH)

    # predicted values
    pv_df = pd.read_csv('data/output/doric_set_no_model_orca.csv')
    # ground truths
    gt_df = pd.read_csv('data/input/experimental_dataset.csv')

    gt_df['middles'] = gt_df.apply(lambda x: Peak.calc_middle(x.begin, x.end, x.seq_len), axis=1)

    # pv_df.shape[0] >>> gt_df.shape[0]
    pv_df = pv_df[pv_df['accession'].isin(gt_df['accession'])]
    pv_df.reset_index(inplace=True)

    compare(pv_df, gt_df, RFC_model)


def DoriC_full(path_DoriC: str, path_exp_set: str) -> None:
    DoriC = pd.read_csv(path_DoriC)
    exp_set = pd.read_csv(path_exp_set)

    # only one version per accession in dataset, no need to check version
    DoriC_E = DoriC[~DoriC['accession'].isin(exp_set['accession'])]

    X_DoriC_E, y_DoriC_E = load_data_labels_from_df(DoriC_E)

    RFC_params = {
        'n_estimators': 600,
        'min_samples_split': 2,
        'min_samples_leaf': 2,
        'max_features': 'log2',
        'max_depth': 70,
        'bootstrap': False
    }

    print("Training on DoriC_E")
    RFC_model = RandomForestClassifier(random_state=RANDOM_STATE, **RFC_params)
    RFC_model.fit(X_DoriC_E, y_DoriC_E)
    joblib.dump(RFC_model, MODEL_PATH)


def run_DoriC_split() -> None:

    DoriC = pd.read_csv(DATA_PATH_DORIC)
    exp_set = pd.read_csv(DATA_PATH_EXP_SET)

    # only one version per accession in dataset, no need to check version
    DoriC_E = DoriC[~DoriC['accession'].isin(exp_set['accession'])]

    X_DoriC_E = DoriC_E[['accession', 'Z_score', 'G_score', 'D_score', 'total_pot']]
    y_DoriC_E = DoriC_E['correct']

    Full_info_DoriC = load_data("data/input/DoriC_chromosome_circular.csv",
                                "data/input/DoriC_complete_circular.csv").drop(columns=['oric_seq', 'lineage'])
    Full_seq_info = pd.read_csv("data/output/doric_set_no_model_orca.csv")

    partial_process = partial(
        process,
        X_DoriC_E=X_DoriC_E,
        y_DoriC_E=y_DoriC_E,
        Full_info_DoriC=Full_info_DoriC,
        Full_seq_info=Full_seq_info,
        exp_set=exp_set
    )

    # for i in range(101):
    #     partial_process(i)
    thresholds = [i for i in range(101)]
    with mp.Pool(os.cpu_count() - 2) as pool:
        pool.map(partial_process, thresholds)

    fieldnames = [
        "threshold",
        "precision_30",
        "precision_30_std",
        "recall_30",
        "recall_30_std",
        "precision_exp",
        "precision_exp_std",
        "recall_exp",
        "recall_exp_std"
    ]
    bf.merge_csvs(file_folder="data/output/precision_recall",
                  merged_csv="merged_precision_recall.csv", fieldnames=fieldnames, headers=True)


def process(i, X_DoriC_E, y_DoriC_E, Full_info_DoriC, Full_seq_info, exp_set):
    if os.path.exists(f"data/output/precision_recall/threshold_{i}.csv"):
        return
    threshold = i/100
    print(f"Starting with threshold: {threshold}", flush=True)

    skf = StratifiedKFold(n_splits=5, random_state=RANDOM_STATE, shuffle=True)
    list_p_70, list_r_70, list_p_exp, list_r_exp = [], [], [], []

    start = perf_counter()
    for fold_num, (train_index, test_index) in enumerate(skf.split(X_DoriC_E, y_DoriC_E)):
        print(f"Fold {fold_num+1} with threshold {threshold}", flush=True)

        X_train, X_test = X_DoriC_E.iloc[train_index, :], X_DoriC_E.iloc[test_index, :]
        y_train, y_test = y_DoriC_E.iloc[train_index], y_DoriC_E.iloc[test_index]

        # print('training model...')
        model = train_DoriC_split(X_train, y_train)

        # print('calculating...')
        p_70, r_70 = compare_70(pv_df=Full_seq_info[Full_seq_info['accession'].isin(X_test["accession"]) & ~Full_seq_info['accession'].isin(
            exp_set['accession'])], gt_df=Full_info_DoriC, model=model, threshold=threshold)
        # print(f"precision 70 %: {p_70:.5f}, recall 70 %: {r_70:.5f}")
        p_exp, r_exp = compare_70(pv_df=Full_seq_info[Full_seq_info['accession'].isin(
            exp_set['accession'])], gt_df=Full_info_DoriC, model=model, threshold=threshold)
        # print(f"precision exp : {p_exp:.5f}, recall exp : {r_exp:.5f}")
        # print(f"Time for fold : {perf_counter() - start:.3f} s")
        # print()
        start = perf_counter()

        list_p_70.append(p_70)
        list_r_70.append(r_70)
        list_p_exp.append(p_exp)
        list_r_exp.append(r_exp)

    print(flush=True)
    print(f"Final Results with threshold: {threshold}!", flush=True)
    print(f"avg precision 70 % {np.mean(list_p_70):.5f} +/- {np.std(list_p_70):.5f}", flush=True)
    print(f"avg recall 70 %    {np.mean(list_r_70):.5f} +/- {np.std(list_r_70):.5f}", flush=True)
    print(f"avg precision exp  {np.mean(list_p_exp):.5f} +/- {np.std(list_p_exp):.5f}", flush=True)
    print(f"avg recall exp     {np.mean(list_r_exp):.5f} +/- {np.std(list_r_exp):.5f}", flush=True)
    print(flush=True)

    pd.DataFrame({
        "threshold": [threshold],
        "precision_30": [np.mean(list_p_70)],
        "precision_30_std": [np.std(list_p_70)],
        "recall_30": [np.mean(list_r_70)],
        "recall_30_std": [np.std(list_r_70)],
        "precision_exp": [np.mean(list_p_exp)],
        "precision_exp_std": [np.std(list_p_exp)],
        "recall_exp": [np.mean(list_r_exp)],
        "recall_exp_std": [np.std(list_r_exp)]
    }).to_csv(f"data/output/precision_recall/threshold_{i}.csv")


def train_DoriC_split(X_train: pd.DataFrame, y_train: pd.Series) -> RandomForestClassifier:

    RFC_params = {
        'n_estimators': 600,
        'min_samples_split': 2,
        'min_samples_leaf': 2,
        'max_features': 'log2',
        'max_depth': 70,
        'bootstrap': False
    }

    X_train = X_train.drop(columns=['accession'])

    RFC_model = RandomForestClassifier(random_state=RANDOM_STATE, **RFC_params)
    RFC_model.fit(X_train, y_train)
    # joblib.dump(RFC_model, 'data/output/machine_learning/DoriC_70_model.pkl')
    return RFC_model


def compare_70(pv_df: pd.DataFrame, gt_df: pd.DataFrame, model: RandomForestClassifier, threshold: float) -> None:
    """
    Compare the test_set (30% of DoriC) by ORCA to the ground truth oriCs for the corresponding accession.
    If the predicted origin is further than 2.5% of the total genome length removed from the ground truth,
    it gets labelled as a `False` prediction. If it is within this distance, it gets labelled as a `True` prediction.
    """

    gt_df['accession'] = gt_df['nc'].str.split(".").str[0]
    # pv_df['nc'] = pv_df['accession'] + '.' + pv_df['version'].map(str)

    # pv_df.shape[0] <<< gt_df.shape[0]
    gt_df = gt_df[gt_df['accession'].isin(pv_df['accession'])]

    gt_df = gt_df.merge(pv_df, on='accession', how='right')[gt_df.columns.to_list() + ['seq_len']]
    gt_df['middles'] = gt_df.apply(lambda x: Peak.calc_middle(x.oric_start, x.oric_end, x.seq_len), axis=1)

    pv_df.reset_index(inplace=True)
    return compare(pv_df, gt_df, model, threshold=threshold)


def old_compare(pv_df: pd.DataFrame, gt_df: pd.DataFrame, model: RandomForestClassifier) -> tuple[float, float]:
    """
    Compare all possible found oriCs by ORCA to the ground truth oriCs for the corresponding accession.
    If the predicted origin is further than 2.5% of the total genome length removed from the ground truth,
    it gets labelled as a `False` prediction. If it is within this distance, it gets labelled as a `True` prediction.
    """

    print("|", "-"*100, "|")
    progress_percent = 0
    print("| ", end="")

    TP, TN, FP, FN = 0, 0, 0, 0
    for i, sample_original in pv_df.iterrows():
        sample = sample_original.copy()

        sample.dropna(inplace=True)
        total_pot_oriCs = len([x for x in sample.axes[0] if 'oriC_middle_' in x])

        predictions = []

        for j in range(total_pot_oriCs):
            features = sample[['Z_score_' + str(j), 'G_score_' + str(j), 'D_score_' +
                               str(j)]].to_list() + [total_pot_oriCs]
            chance_to_be_correct = model.predict_proba(np.asarray(features).reshape(1, -1)).tolist()[0][1]
            predictions.append(chance_to_be_correct)

        idx_val, max_val = max(enumerate(predictions), key=lambda x: x[1])
        oriC_middle = sample['oriC_middle_' + str(idx_val)]

        if max_val > 0.5:  # ORCA thought this was the most likely true origin
            # for-loop because the ground truth could be a bipartite oriC. We count ORCA's oriC as correct if it found one (or both) or the bipartite oriC sites.
            correctness = False
            for _, gt_sample in gt_df[gt_df['accession'] == sample['accession']].iterrows():
                if Peak.calc_dist_points(oriC_middle, gt_sample['middles'], sample['seq_len']) < sample.seq_len * 0.025:
                    correctness = True
            if correctness:  # ORCA said this was a true oriC and it was right
                TP += 1
            else:  # ORCA said this was a true oriC and it was wrong
                FP += 1
        else:  # ORCA did not find any site in this genome that was deemed a true origin
            FN += 1
            # Always results in a False Negative since every genome has an origin.

        if (pv_df.shape[0] < 101):
            if (i/pv_df.shape[0]) * 100 >= (pv_df.shape[0] / 100) * progress_percent:
                new = int(f'{i/pv_df.shape[0] * 100 + 1:.0f}')
                print("-" * (new - progress_percent), end="")
                progress_percent = new
        else:
            if i == (pv_df.shape[0] // 100) * progress_percent:
                progress_percent += 1
                print("-", end="")
    print(" |")

    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    return precision, recall


def compare(pv_df: pd.DataFrame, gt_df: pd.DataFrame, model: RandomForestClassifier, threshold: float = 0.5):

    def get_best_pot_oriC_idx(mat: np.ndarray) -> int:
        """
        Get the best performing potential oriC based on the model's prediction as well as a tie-breaking method.

        Initially the potential oriC is chosen in which the model has the most "confidence", i.e. which one does
        the model think is the most likely to be the True origin.

        In the case of a tie, we look at the highest score of the most important feature (Z-score) between the
        tied potential oriCs. If this is also results in some form of tie, we check the next one, etc.

        `mat`: 2D array with the following columns: prediction_probability, Z-score, G-score, D-score, total_pot_oriCs
        for each potential oriC found in the genome.
        """

        # swap columns around to put them in order of feature importance:
        # [prediction, Z-score, G-score, D-score, total_pot_oriCs] -> [prediction, Z-score, D-score, G-score, total_pot_oriCs]
        mat[:, 2], mat[:, 3] = mat[:, 3].copy(), mat[:, 2].copy()

        max_idx_options = [j for j in range(mat.shape[0])]
        for i in range(mat.shape[1]):
            max_idx_in_curr_col = np.where(mat[max_idx_options, i] == np.max(mat[max_idx_options, i]))[0].tolist()
            max_idx_options = [max_idx_options[j] for j in max_idx_in_curr_col]
            if len(max_idx_options) == 1:
                return max_idx_options[0]
        return max_idx_options[0]

    def calc_classification(sample_original: pd.Series, threshold: float):
        sample = sample_original.copy()
        sample.dropna(inplace=True)

        total_pot_oriCs = len([x for x in sample.axes[0] if 'oriC_middle_' in x])

        predictions = []

        cols = ['Z_score_' + str(j) for j in range(total_pot_oriCs)] \
            + ['G_score_' + str(j) for j in range(total_pot_oriCs)] \
            + ['D_score_' + str(j) for j in range(total_pot_oriCs)]

        features_no_total = sample[cols].to_numpy().reshape((3, total_pot_oriCs)).T
        features = np.concatenate((features_no_total, np.full(
            shape=(total_pot_oriCs, 1), fill_value=total_pot_oriCs)), axis=1)
        predictions = np.vstack(model.predict_proba(features)[:, 1])  # [probability its False , probability its True]
        p_f = np.concatenate((predictions, features), axis=1)
        best_oriC_idx = get_best_pot_oriC_idx(p_f)
        oriC_middle = sample['oriC_middle_' + str(best_oriC_idx)]

        if p_f[best_oriC_idx, 0] >= threshold:  # ORCA thought this was the most likely true origin
            # for-loop because the ground truth could be a bipartite oriC. We count ORCA's oriC as correct if it found one (or both) or the bipartite oriC sites.
            correctness = False
            for gt_middle in sample['middles']:
                if Peak.calc_dist_points(oriC_middle, gt_middle, sample['seq_len']) < sample.seq_len * 0.025:
                    correctness = True
            if correctness:  # ORCA said this was a true oriC and it was right
                return "TP"
            else:  # ORCA said this was a true oriC and it was wrong
                return "FP"
        else:  # ORCA did not find any site in this genome that was deemed a true origin
            return "FN"
            # Always results in a False Negative since every genome has an origin.

    collapsed_df = gt_df.groupby("accession")["middles"].apply(tuple).reset_index()
    merged_df = collapsed_df.merge(pv_df, on="accession").reset_index()

    partial_calc = partial(calc_classification, threshold=threshold)
    c = merged_df.apply(partial_calc, axis=1).value_counts(sort=False)

    TP = c["TP"] if "TP" in c.keys() else 0
    FP = c["FP"] if "FP" in c.keys() else 0
    FN = c["FN"] if "FN" in c.keys() else 0

    precision = TP / (TP + FP) if TP + FP != 0 else 0
    recall = TP / (TP + FN) if TP + FN != 0 else 0
    return precision, recall


def main() -> None:
    run_DoriC_split()


if __name__ == '__main__':
    main()
