# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 07:29:33 2018

@author: Marcus
"""

from toom_cook import *
from karatsuba import *
from compound_multiplication_algorithms import *

import numpy as np


def test_schoolbook_multiply():
    
    # multiply the corresponding f's and g's
    fs = [[1,2], [8, -3, 0, 5], [0], [1,10,-8]]
    gs = [[1,2], [1], [5,5], [-4,0,3]]
    
    # the true products
    answers = [[1,4,4], [8, -3, 0, 5], [0,0], [-4, -40, 35, 30, -24]] 
    
    for i in range(len(fs)):
        assert schoolbook_multiply(fs[i], gs[i]) == answers[i]
        
def test_multiply():
    for poly_length in range(100, 105):
        for algorithm_list in [[], ['2'], ['3'], ['4'], ['2','4'],
                               ['3','2','2'], ['5', '2']]:
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == multiply(f,g, algorithm_list)
    
def test_improved_single_karatsuba():
    for poly_length in range(100, 105):
        for _ in range(5):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == single_karatsuba(f,g)

def test_optimized_toom_3():
    for poly_length in range(100, 110):
        for _ in range(5):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == optimized_toom_3(f,g)

def test_optimized_toom_4():
    for poly_length in range(120, 130):
        for _ in range(5):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == optimized_toom_4(f,g)

def test_optimized_toom_5():
    for poly_length in range(220, 240):
        for _ in range(3):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == optimized_toom_5(f,g)
            
def test_optimized_toom_6():
    for poly_length in range(220, 240):
        for _ in range(3):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == optimized_toom_6(f,g)

def test_karatsuba():
    
    for poly_length in range(100, 110):
        for _ in range(1):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == karatsuba_multiply(f,g)

def test_toom_2():
    
    for poly_length in range(100, 110):
        for _ in range(1):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == toom_cook_multiply(f,g, n=2)
        

def test_toom_3():
    
    for poly_length in range(100, 110):
        for _ in range(1):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == toom_cook_multiply(f,g)
            
            
def test_toom_4():
    
    for poly_length in range(100, 110):
        for _ in range(1):
            f = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            g = [int(x) for x in np.random.randint(-100, 100, size=poly_length)]
            assert schoolbook_multiply(f,g) == toom_cook_multiply(f,g, n=4)