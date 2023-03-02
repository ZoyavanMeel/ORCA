import sys, os, pickle, io

import unittest as ut
import unittest.mock as m

sys.path.insert(0,'..')
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from Handlers import *


class TestFileHandler(ut.TestCase):
    def setUp(self):

        self.acc = "NC_000913.3"
        self.email = "Some@email.address"
        self.out = "Mock/folder"

        with open("NC_000913_3.pkl","rb") as fh:
            seq_rec = pickle.load(fh)
        self.seq_rec = seq_rec
        self.genes = ['dnaA', 'dnaN']

        dnaA = Peak.from_edges(3882325, 3883729, 4641652, 0)
        dnaN = Peak.from_edges(3881220, 3882321, 4641652, 0)
        oriC = Peak.from_edges(3925743, 3925975, 4641652, 0)

        self.seq_dict = {
            'accession': 'NC_000913',
            'version': 3,
            'seq': seq_rec.seq,
            'gene_locations': [
                ('dnaN', dnaN),
                ('dnaA', dnaA)
            ],
            'NCBI_oriC': [('oriC', oriC)],
            'seq_len': 4641652
        }

        self.text_io = io.StringIO("""LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022
                DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
                ACCESSION   NC_000913
                VERSION     NC_000913.3

                gene            190..255
                                /gene="thrL"
                                /locus_tag="b0001"
                                /gene_synonym="ECK0001"
                                /db_xref="ASAP:ABE-0000006"
                                /db_xref="ECOCYC:EG11277"
                                /db_xref="GeneID:944742"
                CDS             190..255
                                /gene="thrL"
                                /locus_tag="b0001"
                                /gene_synonym="ECK0001"
                                /codon_start=1
                                /transl_table=11
                                /product="thr operon leader peptide"
                                /protein_id="NP_414542.1"
                                /db_xref="UniProtKB/Swiss-Prot:P0AD86"
                                /db_xref="ASAP:ABE-0000006"
                                /db_xref="ECOCYC:EG11277"
                                /db_xref="GeneID:944742"
                                /translation="MKRISTTITTTITITTGNGAG"
                """)


    def test_fetch_file_1_good(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.return_value = self.text_io
            x = FileHandler.fetch_file(self.acc, self.email, None, "gbwithparts")
            self.assertEqual(self.text_io, x)


    def test_fetch_file_2_bad_request_with_api_key(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=&retmode=text&tool=biopython&email=Some%40email.address"
            code = 400; msg = "Bad Request"; hdrs = "Some headers"; fp = ""
            mock_efetch.side_effect = HTTPError(url, code, msg, hdrs, fp)

            api_message_1 = f' with API_key \'123\''
            api_message_2 = ' and if your API_key is correctly typed and linked to the given email'
            message = f'Unable to fetch accession: \'ACCESSION\' using \'Some@email.address\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.'

            with self.assertRaisesRegex(ValueError, message):
                FileHandler.fetch_file("ACCESSION", self.email, "123", "gbwithparts")


    def test_fetch_file_3_bad_request_no_api_key(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=&retmode=text&tool=biopython&email=Some%40email.address"
            code = 400; msg = "Bad Request"; hdrs = "Some headers"; fp = ""
            mock_efetch.side_effect = HTTPError(url, code, msg, hdrs, fp)

            message = f'Unable to fetch accession: \'ACCESSION\' using \'Some@email.address\'. Please check if the accession is of an existing chromosomal sequence.'

            with self.assertRaisesRegex(ValueError, message):
                FileHandler.fetch_file("ACCESSION", self.email, None, "gbwithparts")


    def test_fetch_file_4_no_connetion(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.side_effect = URLError("My reason")
            message = 'You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.'

            with self.assertRaisesRegex(ConnectionError, message):
                FileHandler.fetch_file("ACCESSION", self.email, "123", "gbwithparts")


    def test_parse_SeqRecord_1_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, self.genes)
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_2_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, ['dnaA'])
        self.seq_dict['gene_locations'].pop(0)
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_3_good(self):
        res = FileHandler.parse_SeqRecord(self.seq_rec, [])
        self.seq_dict['gene_locations'] = []
        self.assertEqual(res, self.seq_dict)


    def test_get_accession_from_gbk(self):
        acc, version = FileHandler.get_accession_from_gbk(self.text_io)
        self.assertEqual("NC_000913", acc)
        self.assertEqual("3", version)

        first_after_call = self.text_io.readline()
        self.assertEqual("LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022\n", first_after_call)


    def test_save_gbk_1_good(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.listdir') as mock_listdir:

            # efetch return value gets closed by save_gbk, so get the value before the call
            to_write = self.text_io.getvalue()
            mock_efetch.return_value = self.text_io
            mock_listdir.return_value = ["NC_000000_1.gbk", "NC_000001_1.gbk", "NC_000002_1.gbk"]

            # Call the function
            FileHandler.save_gbk(self.acc, self.email, self.out)

            # Assert that the file was opened and written to correctly
            path = os.path.join(self.out, 'NC_000913_3.gbk')
            mock_open.assert_called_once_with(path, 'w')
            mock_open().write.assert_called_once_with(to_write)


    def test_save_gbk_2_file_already_exists(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.listdir') as mock_listdir:

            mock_efetch.return_value = self.text_io
            mock_listdir.return_value = ["NC_000913_3.gbk", "NC_000001_1.gbk", "NC_000002_1.gbk"]

            # Call the function
            message = f'\'NC_000913_3.gbk\' already exists in: {self.out}'
            with self.assertRaisesRegex(FileExistsError, message):
                FileHandler.save_gbk(self.acc, self.email, self.out)


    def test_save_pkl_1_good(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.listdir') as mock_listdir, \
            m.patch('Bio.SeqIO.read') as mock_seqIO_read, \
            m.patch('pickle.dump') as mock_pickle_dump:

            mock_efetch.return_value = self.text_io
            mock_listdir.return_value = ["NC_000000_1.pkl", "NC_000001_1.pkl", "NC_000002_1.pkl"]
            mock_seqIO_read.return_value = self.seq_rec

            # Call the function
            FileHandler.save_pkl(self.acc, self.email, self.out)

            # Assert that the file was opened and written to correctly
            path = os.path.join(self.out, 'NC_000913_3.pkl')
            mock_open.assert_called_once_with(path, 'w')
            mock_pickle_dump.assert_called_once_with(self.seq_rec, mock_open.return_value)


    def test_save_pkl_2_file_already_exists(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.listdir') as mock_listdir, \
            m.patch('Bio.SeqIO.read') as mock_seqIO_read, \
            m.patch('pickle.dump') as mock_pickle_dump:
    
            mock_efetch.return_value = self.text_io
            mock_listdir.return_value = ["NC_000913_3.pkl", "NC_000001_1.pkl", "NC_000002_1.pkl"]

            # Call the function
            message = f'\'NC_000913_3.pkl\' already exists in: {self.out}'
            with self.assertRaisesRegex(FileExistsError, message):
                FileHandler.save_pkl(self.acc, self.email, self.out)


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