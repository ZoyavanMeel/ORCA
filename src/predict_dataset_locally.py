# Standard imports
import multiprocessing as mp
from functools import partial
import os, sys, csv, joblib
import numpy as np

# Set these before running
DATASET        = ['NC_002947', 'NC_002971']
CSV_OUT_FOLDER = 'exp_csvs'
PARALLEL       = False
CPUS           = os.cpu_count() - 1
MAX_ORICS      = 10 # Assumption: No more than 10 oriC for a single organism are predicted

EMAIL   = 'no_need_for_a_real@email_address.com'
API_KEY = None
MODEL   = joblib.load('Machine_Learning/exp_train_model.pkl')

# # Cluster path
# sys.path.append('../OriC_Finder/')

from ORCA_1 import find_oriCs

def prep_prediction(accession, email, api_key, model, csv_path, max_oriCs):

    # Quick check to see if th FASTA has already been processed
    if os.path.exists(csv_path + '/' + accession + '.csv'):
        return

    # preferred_properties, all_oriCs = find_oriCs(sample_path)
    preferred_properties = find_oriCs(accession=accession, email=email, api_key=api_key, model=model)

    row = []
    RefSeq = preferred_properties['accession'].split('.')[0] # RefSeq = RefSeq Accession Number, removing Version Number
    row.append(RefSeq) 
    row.append(preferred_properties['seq_len'])

    # Quality of Prediction processing
    row.append(preferred_properties['gc_conc'])

    # OriC processing
    for i in range(max_oriCs):
        row.append(preferred_properties['oriC_middles'][i]) if i < len(preferred_properties['oriC_middles']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['Z_scores'][i]) if i < len(preferred_properties['Z_scores']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['G_scores'][i]) if i < len(preferred_properties['G_scores']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['D_scores'][i]) if i < len(preferred_properties['D_scores']) else row.append(np.nan)
    for i in range(max_oriCs):
        row.append(preferred_properties['Predictions'][i]) if i < len(preferred_properties['Predictions']) else row.append(np.nan)

    with open(csv_path + '/' + RefSeq + '.csv', 'w') as fh:
        writer = csv.writer(fh)
        writer.writerow(row)
        fh.close()

if __name__ == '__main__':

    if not PARALLEL:
        for accession in DATASET:
            prep_prediction(accession, EMAIL, API_KEY, MODEL, CSV_OUT_FOLDER, MAX_ORICS)
    else:
        prepped_prediction = partial(prep_prediction, email=EMAIL, api_key=API_KEY, csv_path=CSV_OUT_FOLDER, max_oriCs=MAX_ORICS, model=MODEL)
        with mp.Pool(CPUS) as pool:
            pool.map(prepped_prediction, DATASET)