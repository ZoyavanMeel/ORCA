import os, io, pickle

import unittest as ut
import unittest.mock as m

from context import *


class TestCurveProcessing(ut.TestCase):
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
        res_x_min = sorted(CurveProcessing.process_curve(self.x, 'min', 2))
        res_x_max = sorted(CurveProcessing.process_curve(self.x, 'max', 2))
        res_y_min = sorted(CurveProcessing.process_curve(self.y, 'min', 2))
        res_y_max = sorted(CurveProcessing.process_curve(self.y, 'max', 2))
        res_gc_min = sorted(CurveProcessing.process_curve(self.gc, 'min', 2))
        res_gc_max = sorted(CurveProcessing.process_curve(self.gc, 'max', 2))

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