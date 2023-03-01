import sys, os, pickle

import unittest as ut
from unittest.mock import patch

from Bio.SeqFeature import SimpleLocation

sys.path.insert(0,'..')
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from Handlers import *


class TestFileHandler(ut.TestCase):
    def setUp(self):
        with open("NC_000913.pkl","rb") as fh:
            seq_rec = pickle.load(fh)
        self.seq_rec = seq_rec
        self.genes = ['dnaA', 'dnaN']

        dnaN = SimpleLocation(3881220, 3882321, -1)
        dnaA = SimpleLocation(3882325, 3883729, -1)
        oriC = SimpleLocation(3925743, 3925975, +1)

        self.seq_dict = {
            'accession': 'NC_000913',
            'version': 3,
            'seq': seq_rec.seq,
            'dnaA': dnaA,
            'dnaN': dnaN,
            'NCBI_oriC': oriC
        }


    def test_fetch_file_1_good(self):
        with patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.return_value = "TextIO Object"
            x = FileHandler.fetch_file("NC_000913", "Some@email.address", None, "gbwithparts")
            self.assertEqual("TextIO Object", x)


    def test_fetch_file_2_bad_request(self):
        with patch('Bio.Entrez.efetch') as mock_efetch:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=&retmode=text&tool=biopython&email=Some%40email.address"
            code = 400; msg = "Bad Request"; hdrs = "Some headers"; fp = ""
            mock_efetch.side_effect = HTTPError(url, code, msg, hdrs, fp)

            api_message_1 = f' with API_key \'123\''
            api_message_2 = ' and if your API_key is correctly typed and linked to the given email'
            message = f'Unable to fetch accession: \'ACCESSION\' using \'Some@email.address\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.'

            with self.assertRaisesRegex(ValueError, message):
                FileHandler.fetch_file("ACCESSION", "Some@email.address", "123", "gbwithparts")


    def test_fetch_file_3_no_connetion(self):
        with patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.side_effect = URLError("My reason")
            message = 'You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.'

            with self.assertRaisesRegex(ConnectionError, message):
                FileHandler.fetch_file("ACCESSION", "Some@email.address", "123", "gbwithparts")


    def test_parse_SeqRecord_1_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, self.genes)
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_2_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, ['dnaA'])
        self.seq_dict.pop('dnaN')
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_3_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, [])
        self.seq_dict.pop('dnaN')
        self.seq_dict.pop('dnaA')
        self.assertEqual(res, self.seq_dict)



class TestSequenceHandler(ut.TestCase):
    def test_seq_function_1(self):
        # Test case code here
        self.assertEqual(2, 2)

    def test_seq_function_2(self):
        # Test case code here
        self.assertEqual(2, 2)


class TestGeneHandler(ut.TestCase):
    def test_gene_function_1(self):
        # Test case code here
        self.assertEqual(2, 2)


if __name__ == '__main__':
    ut.main()