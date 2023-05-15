# Standard imports
import multiprocessing as mp
from functools import partial
import os, sys, csv, joblib
import numpy as np
from typing import Union

# Set these before running
INPUT_PATH     = 'data\input\dataset.csv'
CSV_OUT_FOLDER = 'data\output\experimental_set_predictions_orca'
OUT_FILE_PATH  = "data\output\experimental_set_predictions_orca\orca.csv"
PARALLEL       = True
CPUS           = os.cpu_count() - 1
MAX_ORICS      = 10 # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'no_need_for_a_real@email_address.com'
API_KEY = None
MODEL   = joblib.load('Machine_Learning/exp_train_model.pkl')

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from ORCA import ORCA
import BioFile

def predict_accession_to_csv(accession: str, email: str, api_key: Union[str, None], model, csv_path: str, max_oriCs: int):
    '''Predict a single `accession` and output the result in a single line CSV file.'''
    # Quick check to see if the accession has already been processed
    out_file = csv_path + '/' + accession.split(".")[0] + '_' + accession.split(".")[1] + '.csv'
    if os.path.exists(out_file):
        return

    orca = ORCA.from_accession(accession, email=email, api_key=api_key, model=model)
    orca.find_oriCs(False, False)

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
    for i in range(max_oriCs):
        row.append(orca.predictions[i]) if i < len(orca.predictions) else row.append(np.nan)

    with open(out_file, 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()


def predict_dataset(input_path: str,
            csv_out_folder: str,
            out_file_path: str,
            parallel: bool,
            cpus: Union[int, None],
            max_oriCs: int,
            email: str,
            api_key: Union[str, None],
            model):
    accessions = []
    with open(input_path, 'r') as fh:
        for row in csv.DictReader(fh):
            accessions.append(row['accession'] + '.' + row['version'])
    accessions = set(accessions)

    if parallel:
        prepped_prediction = partial(predict_accession_to_csv, email=email, api_key=api_key, model=model, csv_path=csv_out_folder, max_oriCs=max_oriCs)
        with mp.Pool(cpus) as pool:
            pool.map(prepped_prediction, accessions)
    else:
        for accession in accessions:
            predict_accession_to_csv(accession, email, api_key, model, csv_out_folder, max_oriCs)

    fieldnames = ['accession', 'version', 'seq_len',]
    fieldnames.extend([f'oriC_middle_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'Z_score_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'G_score_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'D_score_{i}' for i in range(max_oriCs)])
    fieldnames.extend([f'prediction_{i}' for i in range(max_oriCs)])
    BioFile.merge_csvs(csv_out_folder, out_file_path, fieldnames)


if __name__ == '__main__':
    predict_dataset(INPUT_PATH, CSV_OUT_FOLDER, OUT_FILE_PATH, PARALLEL, CPUS, MAX_ORICS, EMAIL, API_KEY, MODEL)