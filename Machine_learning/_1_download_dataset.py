from typing import Iterable
import multiprocessing as mp
from functools import partial
import pandas as pd

import BioFile as bf


INPUT_PATHS = ['data\input\DoriC_chromosome_circular.csv', 'data\input\DoriC_complete_circular.csv']
OUTPUT_PATH = bf.comp_path('data\output\doric_set_no_model_orca_pkl')

EMAIL   = 'zoyavanmeel@gmail.com'
API_KEY = '795d705fb638507c9b2295c89cc64ee88108'
CPUS    = 20


def load_data(*paths) -> pd.DataFrame:
    '''Load data of given paths and concatenates them.'''
    df = pd.DataFrame()
    for path in paths:
        df = pd.concat(objs=[df, pd.read_csv(bf.comp_path(path))], ignore_index=True, join='outer')
    return df


def download_set(data: Iterable[str], output_path: str, cpus: int, email: str, api_key: str = None) -> None:
    downloader = partial(bf.save_pkl, email=email, api_key=api_key, output_folder=output_path)
    with mp.Pool(cpus) as pool:
        pool.map(downloader, data)


if __name__ == "__main__":
    data = set(
        load_data(*INPUT_PATHS)
            .nc # column name
            .astype('string')
            .tolist()
        )
    download_set(data=data, output_path=OUTPUT_PATH, cpus=CPUS, email=EMAIL, api_key=API_KEY)
