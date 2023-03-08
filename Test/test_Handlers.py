import os, io, pickle

import unittest as ut
import unittest.mock as m

from context import *


class TestFileHandler(ut.TestCase):
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
            x = FileHandler.fetch_file(self.acc, self.email, None, "gbwithparts")
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
                FileHandler.fetch_file("ACCESSION", self.email, "123", "gbwithparts")
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
                FileHandler.fetch_file("ACCESSION", self.email, None, "gbwithparts")
            except HTTPError as e:
                self.assertRegex(str(e), message)
                self.assertRegex(e.__notes__[0], note)


    def test_fetch_file_4_no_connetion(self):
        with m.patch('Bio.Entrez.efetch') as mock_efetch:
            mock_efetch.side_effect = URLError("My reason")
            note = 'You are fetching a file from the NCBI servers. Please make sure you have an internet connection to do so.'

            try:
                FileHandler.fetch_file("ACCESSION", self.email, "123", "gbwithparts")
            except URLError as e:
                self.assertRegex(str(e), "My reason")
                self.assertRegex(e.__notes__[0], note)


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
            message = f'\'NC_000913_3.gbk\' already exists in: {self.inp}'
            with self.assertRaisesRegex(FileExistsError, message):
                FileHandler.save_gbk(self.acc, self.email, self.inp)


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
            message = f'\'NC_000913_3.pkl\' already exists in: {self.inp}'
            with self.assertRaisesRegex(FileExistsError, message):
                FileHandler.save_pkl(self.acc, self.email, self.inp)


class TestCurveHandler(ut.TestCase):
    def setUp(self):
        '''
        x = [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 19 18 17
        16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0  1  2  3  4  5  6  7
        8]

        y = [ 0  1  2  3  4  5  6  7  8  9 10  9  8  7  6  5  6  7  8  9 10 11 12 13
        14 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1  0  0  1  2  3  4  5  6
        7  8]

        gc = [  0  -1  -2  -3  -4  -5  -6  -7  -8  -9 -10 -11 -12 -13 -14 -15 -16 -17
        -18 -19 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10  -9  -8  -7  -6  -5
        -6  -7  -8  -9 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1]
        '''

        x = [i for i in range(20)] + [i for i in range(20, 0, -1)] + [i for i in range(0, 10)]
        y = [i for i in range(10)] + [i for i in range(10, 5, -1)] + [i for i in range(5, 15)] + [i for i in range(15, -1, -1)] + [i for i in range(10)]
        gc = [i for i in range(0, -20, -1)] + [i for i in range(-20, -5)] + [i for i in range(-5, -10, -1)] + [i for i in range(-10, 1)]

        self.x = np.asarray(x)
        self.y = np.asarray(y)
        self.gc = np.asarray(gc)


    def test_process_curve_1(self):
        res_x_min = sorted(CurveHandler.process_curve(self.x, 'min', 2))
        res_x_max = sorted(CurveHandler.process_curve(self.x, 'max', 2))
        res_y_min = sorted(CurveHandler.process_curve(self.y, 'min', 2))
        res_y_max = sorted(CurveHandler.process_curve(self.y, 'max', 2))
        res_gc_min = sorted(CurveHandler.process_curve(self.gc, 'min', 2))
        res_gc_max = sorted(CurveHandler.process_curve(self.gc, 'max', 2))

        exp_x_min = [Peak(0, 50, 2), Peak(40, 50, 2)]
        exp_x_max = []
        exp_y_min = []
        exp_y_max = []
        exp_gc_min = []
        exp_gc_max = []

        self.assertEqual(exp_x_min, res_x_min)
        # self.assertEqual(exp_x_max, res_x_max)
        # self.assertEqual(exp_y_min, res_y_min)
        # self.assertEqual(exp_y_max, res_y_max)
        # self.assertEqual(exp_gc_min, res_gc_min)
        # self.assertEqual(exp_gc_max, res_gc_max)
    

    def test_curve_combinations_1(self):
        ...


    def test_detect_peaks_1(self):
        ...
    

    def test_filter_peaks_1(self):
        ...



if __name__ == '__main__':
    ut.main()