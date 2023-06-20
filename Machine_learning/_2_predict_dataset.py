# Standard imports
import multiprocessing as mp
from functools import partial
import os, sys, csv, joblib
import numpy as np
import pandas as pd
from typing import Union

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from ORCA import ORCA
from _1_download_dataset import load_data
import BioFile as bf


# Set these before running
INPUT_PATH     = bf.comp_path('data\output\doric_set_no_model_orca_pkl')
CSV_OUT_FOLDER = bf.comp_path('data\output\doric_set_no_model_orca_csvs')
OUT_FILE_PATH  = bf.comp_path('data\output\doric_set_no_model_orca.csv')
PARALLEL       = True
CPUS           = os.cpu_count() - 1
MAX_ORICS      = 10 # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'no_need_for_a_real@email_address.com'
API_KEY = None
MODEL   = None # joblib.load('Machine_Learning/exp_train_model.pkl')


def predict_pkl_to_csv(path: str, csv_path: str, max_oriCs: int):
    '''Predict a single `accession` and output the result in a single line CSV file.'''

    # Quick check to see if the accession has already been processed
    accession = os.path.split(path)[-1]
    out_file = csv_path + '/' + accession.split("_")[0] + '_' + accession.split("_")[1] + '.csv'
    if os.path.exists(out_file):
        return

    orca = ORCA.from_pkl(path)
    orca.find_oriCs()

    row = []
    row.append(orca.accession) 
    row.append(orca.version)
    row.append(orca.seq_len)

    # OriC processing
    for i in range(max_oriCs):
        row.append(orca.oriC_middles[i]) if i < len(orca.oriC_middles) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(orca.Z_scores[i]) if i < len(orca.Z_scores) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(orca.G_scores[i]) if i < len(orca.G_scores) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(orca.D_scores[i]) if i < len(orca.D_scores) else row.append(np.nan)
    # for i in range(max_oriCs):
    #     row.append(orca.predictions[i]) if i < len(orca.predictions) else row.append(np.nan)

    with open(out_file, 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()


def predict_dataset(
        path_folder: str,
        csv_out_folder: str,
        out_file_path: str,
        parallel: bool,
        cpus: Union[int, None],
        max_oriCs: int
    ):

    paths = os.listdir(path_folder)
    if parallel:
        prepped_prediction = partial(predict_pkl_to_csv, csv_path=csv_out_folder, max_oriCs=max_oriCs)
        with mp.Pool(cpus) as pool:
            pool.map(prepped_prediction, paths)
    else:
        for path in paths:
            predict_pkl_to_csv(path, csv_out_folder, max_oriCs)

    fieldnames = ['accession', 'version', 'seq_len',]
    fieldnames.extend([f'oriC_middle_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'Z_score_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'G_score_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'D_score_{i}' for i in range(max_oriCs)])
    # fieldnames.extend([f'prediction_{i}' for i in range(max_oriCs)])
    bf.merge_csvs(csv_out_folder, out_file_path, fieldnames)


if __name__ == '__main__':
    predict_dataset(
        path_folder=INPUT_PATH,
        csv_out_folder=CSV_OUT_FOLDER,
        out_file_path=OUT_FILE_PATH,
        parallel=PARALLEL, 
        cpus=CPUS,
        max_oriCs=MAX_ORICS
    )