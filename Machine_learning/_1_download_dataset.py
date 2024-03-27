"""Script for downloading all information for the given set of accession numbers and stores them as pickle files of Bio.SeqRecord objects."""

import BioFile as bf
import os
import sys
from typing import Iterable
import multiprocessing as mp
from functools import partial
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))


CLUSTER = "/tudelft.net/staff-umbrella/GeneLocations/ZoyavanMeel/ORCA/"

INPUT_PATHS = [CLUSTER + 'data/input/DoriC_chromosome_circular.csv', CLUSTER + 'data/input/DoriC_complete_circular.csv']
OUTPUT_PATH = CLUSTER + 'data/output/doric_set_no_model_orca_pkl'

EMAIL = 'zoyavanmeel@gmail.com'
API_KEY = 'api_key'
CPUS = 32


def load_data(*paths) -> pd.DataFrame:
    '''Load data of given paths and concatenates them.'''
    df = pd.DataFrame()
    for path in paths:
        df = pd.concat(objs=[df, pd.read_csv(path)], ignore_index=True, join='outer')
    return df


def save_pkl_wrapper(accession: str, email: str, api_key: str, output_folder: str) -> None:
    try:
        bf.save_pkl(accession=accession, email=email, api_key=api_key, output_folder=output_folder)
    except FileExistsError:
        return


def download_set(accessions: Iterable[str], output_path: str, cpus: int, email: str, api_key: str = None) -> None:
    downloader = partial(save_pkl_wrapper, email=email, api_key=api_key, output_folder=output_path)
    with mp.Pool(cpus) as pool:
        pool.map(downloader, accessions)


if __name__ == "__main__":
    accessions = load_data(*INPUT_PATHS)["nc"].astype('string').to_list()
    download_set(accessions=accessions, output_path=OUTPUT_PATH, cpus=CPUS, email=EMAIL, api_key=API_KEY)
