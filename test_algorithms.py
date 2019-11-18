# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 16:39:34 2019

@author: Marcus

This is for testing algorithms.py and nothing else
"""

from algorithms import *
from algorithms_mod import *

import numpy as np

def test_multiply_2048():
    
    for poly_length in range(150, 160):
        for algorithm in [[], [2], [3], [2,2,2,2], [4], 2, 3, [4,2],
                               [3,2,2], [3,3]]:
            f = [int(x) for x in np.random.randint(0, 2048,
                                                    size=poly_length)]
            g = [int(x) for x in np.random.randint(0, 2048, 
                                                    size=poly_length)]
            
            assert schoolbook_mod(f,g, 2048) == multiply_2048(f, g, algorithm)


def test_schoolbook_multiply():
    
    # multiply the corresponding f's and g's
    fs = [[1,2], [8, -3, 0, 5], [0], [1,10,-8]]
    gs = [[1,2], [1], [5,5], [-4,0,3]]
    
    # the true products
    answers = [[1,4,4], [8, -3, 0, 5], [0,0], [-4, -40, 35, 30, -24]] 
    
    for i in range(len(fs)):
        assert schoolbook(fs[i], gs[i]) == answers[i]


def test_algorithm_lists():
    
    for poly_length in range(150, 160):
        for algorithm_list in [[], [2], [3], [4], [2,4], [2,2,2,2],
                               [3,2,2], [5,2], [5], [6], [2,6],
                               [4,5], [8], [8,2], [12]]:
            f = [int(x) for x in np.random.randint(-100, 100,
                                                    size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, 
                                                    size=poly_length)]
            assert schoolbook(f,g) == multiply(f, g, algorithm_list)
            
#def test_single_algorithms():
#    
#    for poly_length in range(150, 160):
#        for algorithm in [2, 3, 4, 5, 6]:
#            for minsize in [2]:#, 6, 10, 50]:
#                f = [int(x) for x in np.random.randint(-100, 100, 
#                                                       size=poly_length)]
#                g = [int(x) for x in np.random.randint(-100, 100, 
#                                                       size=poly_length)]
#                assert schoolbook(f,g) == multiply(f, g, algorithm, 
#                                                   minsize=minsize)



