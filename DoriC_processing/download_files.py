import sys, os
import multiprocessing as mp
from functools import partial

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from old_orca.helper_functions import download

EMAIL    = 'youremail@here.com'
API_KEY  = None # API_key highly recommeneded

ACCESSIONS_FILE = 'DoriC_processing/DoriC_accessions.txt' # accessions that are circular and bacterial
OUTPUT_FOLDER   = 'your_output_folder_here'

# Max. requests = 10 per second, so leave mp.Pool at 10 of you will get blocked.

accessions = []
with open(ACCESSIONS_FILE, 'r') as fh:
    for line in fh:
        accession = line.rstrip('\n')
        accessions.append(accession)

if __name__ == '__main__':
    p_download = partial(download, output_folder=OUTPUT_FOLDER, email=EMAIL, api_key=API_KEY)
    for i in range(len(accessions)):
        download(accessions[i], output_folder=OUTPUT_FOLDER)
        print(f'downloaded: {i+1/len(accessions)}: {accessions[i]}')
    # with mp.Pool(10) as pool:
    #    pool.map(p_download, names)