import unittest as ut

from context import *

# remove later
class TestPeak(ut.TestCase):
    def setUp(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / c    a b           f          h                                                     g    d        e
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        '''
        self.a = Peak(60, 1000, 60)
        self.b = Peak(80, 1000, 60)
        self.c = Peak(10, 1000, 60)
        self.d = Peak(900, 1000, 60)
        self.e = Peak(990, 1000, 60)

        self.f = Peak(200, 1000, 60)
        self.g = Peak(850, 1000, 60)
        self.h = Peak(310, 1000, 60)

        self.lst1 = [self.a, self.b, self.c, self.d, self.e]
        self.lst2 = [self.f, self.g, self.h]


    def test_init_1(self):
        x = Peak(0, 0, 0)
        self.assertIsInstance(x, Peak)


    def test_init_2(self):
        x = Peak(8, 10, 8)
        self.assertIsInstance(x, Peak)


    def test_from_calc_middle(self):
        x = Peak(5, 10, 2)
        y = Peak.from_calc_middle(3, 7, 10, 2)
        self.assertEqual(x, y)


    def test_from_edges(self):
        x = Peak(5, 10, 2)
        y = Peak.from_edges(4, 6, 10)
        self.assertEqual(x, y)


    def test_calc_dist_points_1(self):
        res = Peak.calc_dist_points(self.a.middle, self.b.middle, self.a.seq_len)
        self.assertEqual(res, 20)

    
    def test_calc_dist_points_2(self):
        res = Peak.calc_dist_points(self.a.middle, self.d.middle, self.a.seq_len)
        self.assertEqual(res, 160)


    def test_calc_dist_1(self):
        res = self.a.calc_dist(self.b)
        self.assertEqual(res, 20)

    
    def test_calc_dist_2(self):
        res = self.a.calc_dist(self.d.middle)
        self.assertEqual(res, 160)


    def test_get_peaks_to_merge(self):
        res = Peak.get_peaks_to_merge(self.lst1)
        exp = [(self.a, self.b), (self.c, self.e)]
        self.assertEqual(res, exp)


    def test_match_peaks(self):
        x = [self.a, self.b, self.c]
        y = [self.c, self.d, self.e]

        res = sorted(Peak.match_peaks(x, y))
        exp = sorted([(self.c, self.c), (self.a, self.c), (self.c, self.e)])
        self.assertEqual(res, exp)


    def test_get_adjacency_matrix_1(self):
        res = Peak.get_adjacency_matrix(self.lst1)
        exp = np.array([
            [  0,  20,  50, 160,  70],
            [ 20,   0,  70, 180,  90],
            [ 50,  70,   0, 110,  20], 
            [160, 180, 110,   0,  90],
            [ 70,  90,  20,  90,   0]
        ])
        self.assertTrue(np.array_equal(res, exp))


    def test_get_adjacency_matrix_2(self):
        res = Peak.get_adjacency_matrix(self.lst1, self.lst2)
        exp = np.array([
            #  f,   g,   h
            [140, 210, 250], # a
            [120, 230, 230], # b
            [190, 160, 300], # c
            [300,  50, 410], # d
            [210, 140, 320]  # e
        ])
        self.assertTrue(np.array_equal(res, exp))


    def test_get_adjacency_matrix_3(self):
        lst1 = [x.middle for x in self.lst1]
        lst2 = [x.middle for x in self.lst2]
        res = Peak.get_adjacency_matrix(lst1, lst2, seq_len=1000)
        exp = np.array([
            #  f,   g,   h
            [140, 210, 250], # a
            [120, 230, 230], # b
            [190, 160, 300], # c
            [300,  50, 410], # d
            [210, 140, 320]  # e
        ])
        self.assertTrue(np.array_equal(res, exp))
    

    def test_get_adjacency_matrix_4(self):
        lst2 = [x.middle for x in self.lst2]
        with self.assertRaisesRegex(ValueError, "Elements is `peaks_a` are not the same type as `peaks_b`."):
            Peak.get_adjacency_matrix(self.lst1, lst2)


    def test_get_adjacency_matrix_5(self):
        lst2 = [self.f, self.g.middle, self.h]
        with self.assertRaisesRegex(ValueError, "Not all elements in `peaks_b` are of the same type."):
            Peak.get_adjacency_matrix(self.lst1, lst2)


    def test_get_adjacency_matrix_5(self):
        lst1 = [self.f, self.g.middle, self.h]
        with self.assertRaisesRegex(ValueError, "Not all elements in `peaks_a` are of the same type."):
            Peak.get_adjacency_matrix(lst1)


    def test_get_adjacency_matrix_7(self):
        lst1 = [x.middle for x in self.lst2]
        with self.assertRaisesRegex(ValueError, "Provided list of integers, but did not provide `seq_len`."):
            Peak.get_adjacency_matrix(lst1)


    def test_get_adjacency_matrix_8(self):
        lst2 = [self.f, self.g.middle, self.h]
        with self.assertRaisesRegex(ValueError, "Not all elements in `peaks_b` are of the same type."):
            Peak.get_adjacency_matrix(self.lst1, lst2)


    def test_get_connected_components_idx_1(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / i   i i                                         jj   j     j                                      i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /________                                         ____________                                      __
        '''

        i1 = Peak(  0, 1000, 60)
        i2 = Peak( 40, 1000, 60)
        i3 = Peak( 60, 1000, 60)
        i4 = Peak(990, 1000, 60)

        j1 = Peak(490, 1000, 60)
        j2 = Peak(500, 1000, 60)
        j3 = Peak(540, 1000, 60)
        j4 = Peak(600, 1000, 60)

        peaks = [i1, i2, i3, i4, j1, j2, j3, j4]
        mat = Peak.get_adjacency_matrix(peaks)

        res = sorted(Peak.get_connected_components_idx([i for i in range(len(peaks))], mat, 120))
        exp = sorted([[0, 1, 2, 3], [4, 5, 6, 7]])
        self.assertEqual(res, exp)


    def test_get_connected_components_idx_2(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        /i  i  i  i  i  i  i  i  i  i  i  i  i  i  i  i        i  i  i  i  i  i  i  i  i  i  i  i  i  i  i  i 
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /______________________________________________        _______________________________________________
        '''
        lst_peaks = [Peak(i, 1000, 60) for i in range(0, 1000, 20)]
        indexes = [i for i in range(len(lst_peaks))]

        mat = Peak.get_adjacency_matrix(lst_peaks)

        res = Peak.get_connected_components_idx(indexes, mat, 50)
        exp = [indexes]
        self.assertEqual(res, exp)


    def test_get_connected_components_idx_3(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        /      i  i  i  i  i  i  i  i  i  i  i  i  i  i        i  i  i  i  i  i  i  i  i  i  i  i  i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /      ________________________________________        _____________________________________
        '''
        lst_peaks_1 = [Peak(i, 1000, 60) for i in range(60, 450, 20)]
        lst_peaks_2 = [Peak(i, 1000, 60) for i in range(540, 910, 20)]
        indexes = [i for i in range(len(lst_peaks_1 + lst_peaks_2))]

        mat = Peak.get_adjacency_matrix(lst_peaks_1 + lst_peaks_2)

        res = Peak.get_connected_components_idx(indexes, mat, 50)
        exp = [[i for i in range(len(lst_peaks_1))], [i+len(lst_peaks_1) for i in range(len(lst_peaks_2))]]
        self.assertEqual(res, exp)


    def test_select_connected_groups_2(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / i   i i                                         jj   j     j                                      i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /________                                         ______     _                                      __
        '''
        i1 = Peak(  0, 1000, 60)
        i2 = Peak( 40, 1000, 60)
        i3 = Peak( 60, 1000, 60)
        i4 = Peak(990, 1000, 60)

        j1 = Peak(490, 1000, 60)
        j2 = Peak(500, 1000, 60)
        j3 = Peak(540, 1000, 60)
        j4 = Peak(600, 1000, 60)

        peaks = [i1, i2, i3, i4, j1, j2, j3, j4]
        mat = Peak.get_adjacency_matrix(peaks)
        res = sorted(Peak.select_connected_groups(peaks, mat, 50))
        exp = sorted([[i1, i2, i3, i4], [j1, j2, j3], [j4]])
        self.assertEqual(res, exp)


    def test_select_connected_groups_3(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / i   i i  i   i i                                jj   j   j                                        i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /__   ___  _   ___                                ______   _                                        __
        '''
        i1 = Peak(  0, 1000, 60)
        i2 = Peak( 40, 1000, 60)
        i3 = Peak( 60, 1000, 60)
        i4 = Peak(990, 1000, 60)
        i5 = Peak(100, 1000, 60)
        i6 = Peak(140, 1000, 60)
        i7 = Peak(160, 1000, 60)

        j1 = Peak(490, 1000, 60)
        j2 = Peak(500, 1000, 60)
        j3 = Peak(540, 1000, 60)
        j4 = Peak(600, 1000, 60)

        peaks = [i1, i2, i3, i4, i5, i6, i7, j1, j2, j3, j4]
        mat = Peak.get_adjacency_matrix(peaks)
        # print("-------------------------------------------------")
        res = sorted(Peak.select_connected_groups(peaks, mat, 50))
        # print([[str(i) for i in x] for x in res])
        exp = sorted([[i1, i4], [i2, i3], [i5], [i6, i7], [j1, j2, j3], [j4]])
        # print("-------------------------------------------------")
        self.assertEqual(res, exp)


    def test_select_connected_groups_4(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / i   i i                                jjjjjjjjjj jjjjjjjjjj                                      i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /________                                __________ __________                                      __
        '''
        i1 = Peak(  0, 1000, 60)
        i2 = Peak( 40, 1000, 60)
        i3 = Peak( 60, 1000, 60)
        i4 = Peak(990, 1000, 60)

        jlst1 = [Peak(400+i, 1000, 60) for i in range(0, 100, 10)]
        jlst2 = [Peak(510+i, 1000, 60) for i in range(0, 90, 10)]

        peaks = [i1, i2, i3, i4] + jlst1 + jlst2
        mat = Peak.get_adjacency_matrix(peaks)

        res = sorted(Peak.select_connected_groups(peaks, mat, 50))
        # print([[str(i) for i in x] for x in res])
        exp = sorted([[i1, i2, i3, i4], jlst1, jlst2])
        self.assertEqual(res, exp)


    def test_select_connected_groups_5(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10
        / i   i i                           j j  jjj  jjj j jjj j jj j   j                                  i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /________                           ___  ___  ___ _ ___ _ __ _   _                                  __
        '''
        i1 = Peak(  0, 1000, 60)
        i2 = Peak( 40, 1000, 60)
        i3 = Peak( 60, 1000, 60)
        i4 = Peak(990, 1000, 60)

        js = [350, 370, 400, 410, 420, 450, 460, 470, 490, 510, 520, 530, 550, 570, 580, 600, 640]
        j = [Peak(i, 1000, 60) for i in js]

        peaks = [i1, i2, i3, i4] + j
        mat = Peak.get_adjacency_matrix(peaks)

        res = sorted(Peak.select_connected_groups(peaks, mat, 50))
        # print([[str(i) for i in x] for x in res])
        exp = sorted([[i1, i2, i3, i4], j[0:2], j[2:5], j[5:8], j[8:9], j[9:12], j[12:13], j[13:15], j[15:16], j[16:]])
        # print([[str(i) for i in x] for x in exp])
        self.assertEqual(res, exp)


    def test_select_connected_groups_6(self):
        '''
        /0         1         2         3         4         5         6         7         8         9         10        1         2         3         4         5         6         7         8         9         20
        / i   i i                           j j  jjj  jjj j jjj j jj j   j                                  i i   i i                           j j  jjj  jjj j jjj j jj j   j                                  i
        /|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
        /________                           ___________________________  _                                  _________                           ___________________________  _                                  __
        '''
        i1 = Peak(   0, 2000, 60)
        i2 = Peak(  40, 2000, 60)
        i3 = Peak(  60, 2000, 60)
        i4 = Peak( 990, 2000, 60)

        i5 = Peak(1000, 2000, 60)
        i6 = Peak(1040, 2000, 60)
        i7 = Peak(1060, 2000, 60)
        i8 = Peak(1990, 2000, 60)

        js = [350, 370, 400, 410, 420, 450, 460, 470, 490, 510, 520, 530, 550, 570, 580, 600, 640]
        j1 = [Peak(i, 2000, 60) for i in js]
        j2 = [Peak(i+1000, 2000, 60) for i in js]

        peaks = [i1, i2, i3, i4, i5, i6, i7, i8] + j1 + j2
        mat = Peak.get_adjacency_matrix(peaks)

        # twice as long so 5% = 100
        # would be perfect if all j's were together, but this is a tuning issue, not a testing issue
        # print('======================================================')
        res = sorted(Peak.select_connected_groups(peaks, mat, 100))
        # print([[str(i) for i in x] for x in res])
        # print('======================================================')
        exp = sorted([[i1, i2, i3, i8], [i4, i5, i6, i7], j1[:-1], [j1[-1]], j2[:-1], [j2[-1]]])
        self.assertEqual(res, exp)


    def test_calc_middle_1(self):
        self.assertEqual(Peak.calc_middle(self.a, self.b), 70)


    def test_calc_middle_2(self):
        self.assertEqual(Peak.calc_middle(self.a, self.e), 26)


    def test_calc_middle_3(self):
        self.assertEqual(Peak.calc_middle(self.a, self.g), 955)


    def test_calc_middle_4(self):
        with self.assertRaisesRegex(ValueError, "No curve size found in any parameters. This method is for circular DNA only."):
            Peak.calc_middle(5, 3)


    def test_intersecting_windows_1(self):
        self.assertTrue(self.a.intersecting_windows(self.a))
        self.assertTrue(self.a.intersecting_windows(self.b))
        self.assertTrue(self.a.intersecting_windows(self.c))
        self.assertFalse(self.a.intersecting_windows(self.d))
        self.assertFalse(self.a.intersecting_windows(self.e))
        self.assertFalse(self.b.intersecting_windows(self.e))
        self.assertTrue(self.c.intersecting_windows(self.e))


    def test_intersecting_windows_2(self):
        with self.assertRaisesRegex(ValueError, "other is a <class 'float'>, must be a Peak object"):
            self.a.intersecting_windows(3.5)


    def test_contains_point(self):
        self.assertTrue(self.a.contains_point(40))
        self.assertTrue(self.a.contains_point(self.b))
        self.assertTrue(self.b.contains_point(self.a))
        self.assertFalse(self.a.contains_point(self.d))


    def test_contains(self):
        self.assertTrue(40 in self.a)
        self.assertTrue(self.b in self.a)
        self.assertTrue(self.a in self.b)
        self.assertFalse(self.d in self.a)


    def test_hash(self):
        b = Peak(60, 1000, 60)
        self.assertEqual(hash(self.a), hash(b))
        self.assertNotEqual(hash(self.a), hash(self.b))


    def test_eq(self):
        self.assertEqual(self.a, self.a)
        self.assertEqual(self.a, Peak(60, 1000, 60))
        self.assertNotEqual(self.a, self.b)
        self.assertNotEqual(self.a, 3)


    def test_repr(self):
        s = 'Peak(middle=60, window_size=60, seq_len=1000)'
        self.assertEqual(self.a.__repr__(), s)


    def test_str(self):
        self.assertEqual(str(self.a), "60")


    def test_lt(self):
        self.assertTrue(self.a < self.b)
        self.assertFalse(self.a < self.a)


    def test_gt(self):
        self.assertTrue(self.b > self.a)
        self.assertFalse(self.b > self.b)


    def test_add(self):
        self.assertEqual(Peak(70, 1000, 60), self.a + self.c)
        self.assertEqual(70, self.a + 10)


    def test_sub(self):
        self.assertEqual(Peak(50, 1000, 60), self.a - self.c)
        self.assertEqual(50, self.a - 10)


    def test_radd(self):
        self.assertEqual(160, 100 + self.a)


    def test_rsub(self):
        self.assertEqual(40, 100 - self.a)
    

    def test_neg(self):
        self.assertEqual(Peak(-self.a.middle, self.a.seq_len, self.a.window_size), -self.a)


    def test_float(self):
        self.assertEqual(float(self.a), 60.0)


if __name__ == '__main__':
    ut.main()