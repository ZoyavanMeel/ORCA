import os, io, pickle

import unittest as ut
import unittest.mock as m

from context import *

 
class TestBioFile(ut.TestCase):
    def setUp(self):

        self.acc = "NC_000913.3"
        self.email = "Some@email.address"
        self.inp = "Mock/folder/data/input"
        self.out = "Mock/folder/data/output"

        with open("data/input/NC_000913_3.pkl","rb") as fh:
            seq_rec = pickle.load(fh)
        self.seq_rec = seq_rec
        self.genes = ['dnaA', 'dnaN']

        dnaA = Peak.from_calc_middle(3882325, 3883729, 4641652, 0)
        dnaN = Peak.from_calc_middle(3881220, 3882321, 4641652, 0)
        oriC = Peak.from_calc_middle(3925743, 3925975, 4641652, 0)

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
            x = BioFile.fetch_file(self.acc, self.email, "gbwithparts", None)
            self.assertEqual(self.text_io, x)


    def test_fetch_file_2_bad_request_with_api_key(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=&retmode=text&tool=biopython&email=Some%40email.address"
            code = 400; msg = "Bad Request"; hdrs = "Some headers"; fp = ""
            mock_efetch.side_effect = HTTPError(url, code, msg, hdrs, fp)

            api_message_1 = f' with API_key \'123\''
            api_message_2 = ' and if your API_key is correctly typed and linked to the given email'
            note = f'Unable to fetch accession: \'ACCESSION\' using \'Some@email.address\'{api_message_1}. Please check if the accession is of an existing chromosomal sequence{api_message_2}.'
            message = "HTTP Error 400: Bad Request"

            try:
                BioFile.fetch_file("ACCESSION", self.email, "gbwithparts", "123")
            except HTTPError as e:
                self.assertRegex(str(e), message)
                self.assertRegex(e.__notes__[0], note)


    def test_fetch_file_3_bad_request_no_api_key(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=&retmode=text&tool=biopython&email=Some%40email.address"
            code = 400; msg = "Bad Request"; hdrs = "Some headers"; fp = ""
            mock_efetch.side_effect = HTTPError(url, code, msg, hdrs, fp)

            note = f'Unable to fetch accession: \'ACCESSION\' using \'Some@email.address\'. Please check if the accession is of an existing chromosomal sequence.'
            message = "HTTP Error 400: Bad Request"
    
            try:
                BioFile.fetch_file("ACCESSION", self.email, "gbwithparts", None)
            except HTTPError as e:
                self.assertRegex(str(e), message)
                self.assertRegex(e.__notes__[0], note)


    def test_fetch_file_4_no_connetion(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.side_effect = URLError("My reason")
            note = 'You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.'

            try:
                BioFile.fetch_file("ACCESSION", self.email, "gbwithparts", "123")
            except URLError as e:
                self.assertRegex(str(e), "My reason")
                self.assertRegex(e.__notes__[0], note)


    def test_parse_SeqRecord_1_good(self):
        res = BioFile.parse_SeqRecord(self.seq_rec, self.genes)
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_2_good(self):
        res = BioFile.parse_SeqRecord(self.seq_rec, ['dnaA'])
        self.seq_dict['gene_locations'].pop(0)
        self.assertEqual(res, self.seq_dict)


    def test_parse_SeqRecord_3_good(self):
        res = BioFile.parse_SeqRecord(self.seq_rec, [])
        self.seq_dict['gene_locations'] = []
        self.assertEqual(res, self.seq_dict)


    def test_save_gbk_1_good(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.listdir') as mock_listdir:

            # efetch return value gets closed by save_gbk, so get the value before the call
            to_write_a = """LOCUS       NC_000913            4641652 bp    DNA     circular CON 09-MAR-2022
                DEFINITION  Escherichia coli str. K-12 substr. MG1655, complete genome.
                ACCESSION   NC_000913
                VERSION     NC_000913.3
            """

            to_write_b = """
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
                """
            mock_efetch.return_value = self.text_io
            mock_listdir.return_value = ["NC_000000_1.gbk", "NC_000001_1.gbk", "NC_000002_1.gbk"]

            # Call the function
            BioFile.save_gbk(self.acc, self.email, self.out)

            # Assert that the file was opened and written to correctly
            path = os.path.join(self.out, 'NC_000913_3.gbk')
            mock_open.assert_called_with(path, 'w')
            mock_open().write.has_calls([to_write_a, to_write_b])


    def test_save_gbk_2_file_already_exists(self):
        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.path.exists') as mock_exists:

            mock_efetch.return_value = self.text_io
            mock_exists.return_value = True

            # Call the function
            message = f'\'NC_000913_3.gbk\' already exists in: {self.inp}'
            with self.assertRaisesRegex(FileExistsError, message):
                BioFile.save_gbk(self.acc, self.email, self.inp)


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
            BioFile.save_pkl(self.acc, self.email, self.out)

            # Assert that the file was opened and written to correctly
            path = os.path.join(self.out, 'NC_000913_3.pkl')
            mock_open.assert_called_once_with(path, 'wb')
            mock_pickle_dump.assert_called_once_with(self.seq_rec, mock_open.return_value)


    def test_save_pkl_2_file_already_exists(self):
        class MockSeqIO:
            pass

        # Set up the mock file object
        with m.patch('builtins.open', m.mock_open()) as mock_open, \
            m.patch('Bio.Entrez.efetch') as mock_efetch, \
            m.patch('os.path.exists') as mock_exists, \
            m.patch('Bio.SeqIO.read') as mock_read:
    
            mock_efetch.return_value = self.text_io
            mock_exists.return_value = True
            seqio = MockSeqIO()
            seqio.id = "NC_000913.3"
            mock_read.return_value = seqio

            # Call the function
            message = f'\'NC_000913_3.pkl\' already exists in: {self.inp}'
            with self.assertRaisesRegex(FileExistsError, message):
                BioFile.save_pkl(self.acc, self.email, self.inp)

if __name__ == "__main__":
    ut.main()