import os, shutil, csv
from Bio import SeqIO, Entrez
from typing import TextIO, Union, Tuple
from urllib.error import HTTPError, URLError

class FileHandler:

    @staticmethod
    def fetch_file(accession: str, email: str, api_key: Union[str, None], rettype: str) -> TextIO:
        """Downloads the given file_type of the given accession in temporary memory"""
        Entrez.email = email
        if api_key is not None: Entrez.api_key = api_key
        try:
            return Entrez.efetch(db="nuccore", id=accession, rettype=rettype, retmode="text")
        except HTTPError as BadRequest:
            (api_message_1, api_message_2) = (f' with API_key \'{api_key}\'', ' and if your API_key if correctly typed and linked to the given email') if api_key is not None else ('', '')
            raise ValueError(f'Unable to fetch accession: \'{accession}\' using \'{email}\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.') from BadRequest
        except URLError as NoConnection:
            raise ConnectionError('You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.') from NoConnection


    @staticmethod
    def read_FASTA(handle: Union[TextIO, str]) -> Tuple[str, str]:
        """Read a FASTA file and returns the name and sequence of only the first entry in the file"""
        Seq_records = SeqIO.parse(handle, 'fasta')
        Seq_obj = next(Seq_records)
        return Seq_obj.id, Seq_obj.seq


    @staticmethod
    def move_fastas(db_loc, on_cluster=True, split=4):
        '''
        Split one folder into 'split' amount of folders with roughly the same amount of files in them.
        Used for easier parallel processing. Instead of one big job with 4000 samples. Run 4 jobs with 1000 samples at the same time.
        '''
        if on_cluster:
            path = db_loc + '/bacteria'
            samples = os.listdir(path)
        else:
            path = db_loc + '/chromosomes_only'
            samples = os.listdir(path)

        samples_per_split = len(samples)//split
        for i in range(split):
            new_folder = db_loc + '/' + str(i)
            os.mkdir(new_folder)
            start = 0 + i*samples_per_split
            stop  = (i+1)*samples_per_split if i != split-1 else None
            for sample in samples[start:stop]:
                shutil.move(path + '/' + sample, new_folder + '/' + sample)


    @staticmethod
    def download(accession: str, output_folder: str, email: str = None, api_key: str = None):
        '''Download the proper file types for large dataset analysis into the output_folder.'''
        Entrez.email = email
        if api_key is not None:
            Entrez.api_key = api_key

        if not os.path.exists(output_folder + f'/gene_info_files/{accession}.fasta'):
            check_handle = Entrez.efetch(db='nuccore', id=accession, rettype='gb', retmode='text')
            annotations = next(SeqIO.parse(check_handle, 'genbank').records).annotations
            if annotations['taxonomy'][0].lower() == 'bacteria' and annotations['topology'].lower() == 'circular':
                handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta_cds_na", retmode="text")
                with open(output_folder + f'/gene_info_files/{accession}.fasta', 'w') as fh:
                    lines = handle.read()
                    fh.write(lines)
                    handle.close()
                    fh.close()
                handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
                with open(output_folder + f'/bacteria/{accession}.fasta', 'w') as fh:
                    lines = handle.read()
                    fh.write(lines)
                    handle.close()
                    fh.close()

    @staticmethod
    def merge_csvs(file_folder, merged_csv, fieldnames, length=-1, headers=False):
        '''
        Merge multiple csvs into one csv.
        Arguments:
            file_folder : path to folder with csvs that have to be merged
            merged_csv  : name of the merged csv
            fieldnames  : list of fieldnames for the merged_csv
            length      : amount of rows in each single csv. -1, if the length of each single
                        csv is not known or not the same.
            headers     : if the single csvs have headers or not
        '''
        file_list = os.listdir(file_folder)
        with open(merged_csv, 'w', newline='') as fh_out:
            writer = csv.writer(fh_out)
            writer.writerow(fieldnames)
            for file in file_list:
                with open(file_folder + file, 'r') as fh_in:
                    reader = csv.reader(fh_in)
                    for i, row in enumerate(reader):
                        if headers and i == 0:
                            pass
                        elif i == length:
                            break
                        else:
                            writer.writerow(row)