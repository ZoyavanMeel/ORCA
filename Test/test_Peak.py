import sys, os

import unittest as ut
import unittest.mock as m

sys.path.insert(0,'..')
os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

from Peak import *


class TestPeak(ut.TestCase):
    def setUp(self):
        self.a = Peak(60, 1000, 60)
        self.b = Peak(80, 1000, 60)
        self.c = Peak(10, 1000, 60)
        self.d = Peak(900, 1000, 60)
        self.e = Peak(990, 1000, 60)
        self.lst = [self.a, self.b, self.c, self.d, self.e]

    def test_init(self):
        x = Peak(0, 0, 0)
        self.assertIsInstance(x, Peak)


    def test_from_calc_middle(self):
        x = Peak(5, 10, 2)
        y = Peak.from_calc_middle(3, 7, 10, 2)
        self.assertEqual(x, y)


    def test_from_edges(self):
        x = Peak(5, 10, 2)
        y = Peak.from_edges(4, 6, 10)
        self.assertEqual(x, y)


    def test_calc_dist_1(self):
        res = Peak.calc_dist(self.a.middle, self.b.middle, self.a.seq_len)
        self.assertEqual(res, 20)

    
    def test_calc_dist_2(self):
        res = Peak.calc_dist(self.a.middle, self.d.middle, self.a.seq_len)
        self.assertEqual(res, 160)


    def test_get_peaks_to_merge(self):
        res = Peak.get_peaks_to_merge(self.lst)
        exp = [(self.a, self.b), (self.c, self.e)]
        self.assertEqual(res, exp)


    def test_match_peaks(self):
        x = [self.a, self.b, self.c]
        y = [self.c, self.d, self.e]

        res = sorted(Peak.match_peaks(x, y))
        exp = sorted([(self.c, self.c), (self.a, self.c), (self.c, self.e)])
        self.assertEqual(res, exp)
        


    def test_get_adjacency_matrix(self):
        ...


    def test_get_connected_groups(self):
        ...


    def test_get_middle(self):
        ...


    def test_intersecting_windows(self):
        self.assertTrue(self.a.intersecting_windows(self.a))
        self.assertTrue(self.a.intersecting_windows(self.b))
        self.assertTrue(self.a.intersecting_windows(self.c))
        self.assertFalse(self.a.intersecting_windows(self.d))
        self.assertFalse(self.a.intersecting_windows(self.e))
        self.assertFalse(self.b.intersecting_windows(self.e))
        self.assertTrue(self.c.intersecting_windows(self.e))



    def test_contains_point(self):
        ...


    def test_hash(self):
        ...


    def test_eq(self):
        ...


    def test_repr(self):
        ...


    def test_str(self):
        ...


    def test_lt(self):
        ...


    def test_gt(self):
        ...


    def test_add(self):
        ...


    def test_sub(self):
        ...


    def test_radd(self):
        ...


    def test_rsub(self):
        ...


if __name__ == '__main__':
    ut.main()