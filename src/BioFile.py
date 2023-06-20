"""Module for handing biodata files and whatnot."""

import os, csv, pickle, re

from typing import TextIO, Union
from urllib.error import HTTPError, URLError

from Bio import SeqIO, Entrez

from Peak import Peak


def fetch_file(accession: str, email: str, api_key: Union[str, None], rettype: str) -> TextIO:
    """Downloads the given file_type of the given accession in temporary memory"""
    Entrez.email = email
    if api_key is not None: Entrez.api_key = api_key
    try:
        return Entrez.efetch(db="nuccore", id=accession, rettype=rettype, retmode="text")
    except HTTPError as BadRequest:
        if api_key is not None:
            api_message_1 = f' with API_key \'{api_key}\''
            api_message_2 = ' and if your API_key is correctly typed and linked to the given email'
        else :
            api_message_1, api_message_2 = '', ''
        BadRequest.add_note(f'Unable to fetch accession: \'{accession}\' using \'{email}\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.')
        raise
    except URLError as NoConnection:
        NoConnection.add_note('You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.')
        raise


def parse_SeqRecord(record: SeqIO.SeqRecord, genes_of_interest: list[str]) -> dict:
    """Extracts the sequence and positions of the genes of interest from a SeqRecord."""
    seq_dict = dict()
    accession, version = tuple(record.id.split('.'))
    seq_dict.update({'accession': accession})
    seq_dict.update({'version': int(version)})
    seq_dict.update({'seq': str(record.seq)}) # I only use the sequence to loop over once. strings are much faster for this.
    seq_dict.update({'seq_len': len(record.seq)})
    seq_dict.update({'gene_locations': []})
    seq_dict.update({'NCBI_oriC': []})

    for feature in record.features:
        # is this feature a coding sequence and a gene and is its name something we are looking for?
        if feature.type == 'CDS' and 'gene' in feature.qualifiers and feature.qualifiers['gene'][0] in genes_of_interest:
            gene_loc = Peak.from_calc_middle(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
            seq_dict['gene_locations'].append( (feature.qualifiers['gene'][0], gene_loc) )
        # just in case this SeqRecord has an annotated oriC!
        if feature.type == 'rep_origin':
            oriC = Peak.from_calc_middle(int(feature.location.start), int(feature.location.end), len(record.seq), 0)
            seq_dict['NCBI_oriC'].append( ('oriC', oriC) )
    return seq_dict


def save_gbk(accession: str, email: str, output_folder: str, api_key: str = None):
    '''
    Download the GenBank file of a given accession and save it into the output_folder.
    Name of the file will be: `{accession}_{version}.gbk`.
    Raises a `FileExistsError` if a file with the generated name already exists in the output folder.
    If version is not provided in the accession, then the function downloads the latest version.
    '''
    with fetch_file(accession, email, api_key, rettype="gbwithparts") as fh:
        # Quick search for the version of the accession that was downloaded.
        acc_v = "."
        header = ""
        for line in fh:
            # Underlying stream is not seekable, so cannot reset back to top of file
            # -> need to keep track of what has been read
            header += line
            # Always appears around the top of the file.
            if "VERSION" in line:
                acc_v = [x for x in line.strip("\n").split(" ")][-1]
                break
        acc, version = tuple(acc_v.split("."))

        # Check if a file with the same name already exists
        if os.path.exists(os.path.join(output_folder, acc + '_' + version + '.gbk')):
            fh.close()
            raise FileExistsError(f'\'{acc}_{version}.gbk\' already exists in: {output_folder}')

        # Save contents to path
        file_path = os.path.join(output_folder, acc + '_' + version + '.gbk')
        with open(file_path, 'w') as oh:
            oh.write(header)
            oh.write(fh.read())
            oh.close()
        fh.close()


def save_pkl(accession: str, email: str, output_folder: str, api_key: str = None):
    '''
    Download the GenBank file of a given accession, parses it with Biopython into a SeqRecord, and save it into the output_folder.
    Name of the file will be: `{accession}_{version}.pkl`.
    Raises a `FileExistsError` if a file with the generated name already exists in the output folder.
    If version is not provided in the accession, then the function downloads the latest version.
    '''
    with fetch_file(accession, email, api_key, rettype="gbwithparts") as fh:

        # Parse gbk file
        seq_rec = SeqIO.read(fh, 'gb')
        acc, version = seq_rec.id.split('.')
        
        # Check if a file with the same name already exists
        if os.path.exists(os.path.join(output_folder, acc + '_' + version + '.pkl')):
            fh.close()
            raise FileExistsError(f'\'{acc}_{version}.pkl\' already exists in: {output_folder}')
        
        # Save contents to path
        file_path = os.path.join(output_folder, acc + '_' + version + '.pkl')
        with open(file_path, 'wb') as oh:
            pickle.dump(seq_rec, oh)
            oh.close()
        fh.close()


def merge_csvs(file_folder: str, merged_csv: str, fieldnames: list[str], length: int = -1, headers: bool = False):
    '''
    Merge multiple csvs into one csv.
    Arguments:
    - file_folder : path to folder with csvs that have to be merged
    - merged_csv  : name of the merged csv
    - fieldnames  : list of fieldnames for the merged_csv
    - length      : amount of rows in each single csv. -1, if the length of each single
                    csv is not known or not the same.
    - headers     : if the single csvs have headers or not
    '''
    file_list = os.listdir(file_folder)
    with open(merged_csv, 'w', newline='') as fh_out:
        writer = csv.writer(fh_out)
        writer.writerow(fieldnames)
        for file in file_list:
            with open(os.path.join(file_folder, file), 'r') as fh_in:
                reader = csv.reader(fh_in)
                for i, row in enumerate(reader):
                    if headers and i == 0:
                        pass
                    elif i == length:
                        break
                    else:
                        writer.writerow(row)


def comp_path(path: str) -> str:
    '''Re-join file paths so they are compatible with the current OS.'''
    dirs = re.split(r"\\|\\\\|\/", path)
    if path.startswith(("\\", "\\\\", "\/")):
        return os.sep + os.path.join(*dirs)
    return os.path.join(*dirs)