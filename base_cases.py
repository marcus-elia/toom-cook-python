# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 10:33:29 2019

@author: Marcus

This file is for experimentally determining when schoolbook multiplication
becomes slower than Karatsuba, Toom-3, and Toom-4 on my computer in Python
given my implementations.  Note that my implementation for Karatsuba is fine,
while my implementations for Toom-3 and Toom-4 probably have more matrix 
operations than necessary, and they will be improved if they are bad now.
"""

from karatsuba import *
from toom_cook import *
import numpy as np
import matplotlib.pyplot as plt
import time

def plot_base_case(algorithm, num_trials=500, low_deg=5, high_deg=40, filename=None):
    """ Makes a plot showing when the schoolbook and other curves cross.
        Choose algorithm from 'Single Karatsuba', 'Single Toom-3', or
        'Single Toom-4'"""
    # the degrees we are testing and plotting
    degrees = range(low_deg, high_deg+1)
    
    # dictionaries keeping track of sums of times. We will divide them later
    # to get the averages
    other_times = {degree:0 for degree in degrees}
    sb_times = {degree:0 for degree in degrees}
    
    
    for _ in range(num_trials):
                
        for degree in degrees:
            f = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            g = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            
            # time the schoolbook
            start_time = time.time()
            schoolbook_multiply(f, g)
            sb_times[degree] += time.time() - start_time
            
            # time the other
            if algorithm == 'Single Karatsuba':
                start_time = time.time()
                single_karatsuba(f, g)
                other_times[degree] += time.time() - start_time
                
            elif algorithm == 'Single Toom-3':
                start_time = time.time()
                single_toom_3(f, g)
                other_times[degree] += time.time() - start_time
                
            elif algorithm == 'Single Toom-4':
                start_time = time.time()
                single_toom_4(f, g)
                other_times[degree] += time.time() - start_time
                
            elif algorithm == 'Single Toom-5':
                start_time = time.time()
                single_toom_5(f, g)
                other_times[degree] += time.time() - start_time
                
            elif algorithm == 'Single Toom-6':
                start_time = time.time()
                single_toom_6(f, g)
                other_times[degree] += time.time() - start_time
            
    sb_times = [sb_times[degree]/num_trials for degree in degrees]
    other_times = [other_times[degree]/num_trials for degree in degrees]
        
    plt.plot(degrees, sb_times, lw=3, label='Schoolbook')
    plt.plot(degrees,other_times, lw=3, label=algorithm)
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Time (seconds)", size=15)
    plt.title("Average Time for Polynomial Multiplication", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return

def does_optimization_help(num_trials=100):
    """ Just for myself, to verify that optimizing helped."""
    # the degrees we are testing and plotting
    degrees = range(10, 60)
    
    # dictionaries keeping track of sums of times. We will divide them later
    # to get the averages
    optimized_times = {degree:0 for degree in degrees}
    un_times = {degree:0 for degree in degrees}
    
    
    for _ in range(num_trials):
                
        for degree in degrees:
            f = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            g = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            
            # time the un-optimized
            start_time = time.time()
            toom_cook_multiply(f, g, n=3, minsize=10)
            un_times[degree] += time.time() - start_time
            
            # time the optimized
            start_time = time.time()
            optimized_toom_3(f, g, minsize=10)
            optimized_times[degree] += time.time() - start_time
            
    un_times = [un_times[degree]/num_trials for degree in degrees]
    optimized_times = [optimized_times[degree]/num_trials for degree in degrees]
        
    plt.plot(degrees, un_times, lw=3, label='Not Optimized')
    plt.plot(degrees,optimized_times, lw=3, label='Optimized')
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Time (seconds)", size=15)
    plt.title("Average Time for Polynomial Multiplication", size=15)
    plt.show()
    return

plot_base_case('Single Karatsuba', num_trials=100, low_deg=5, high_deg=50)#, filename='Karatsuba_Base_Case.jpg')
#plot_base_case('Single Toom-3', num_trials=5000, low_deg=10, high_deg=6)#0, filename='Toom_3_Base_Case.jpg')
#plot_base_case('Single Toom-6', num_trials=500, low_deg=60, high_deg=95)#, filename='Toom_5_Base_Case.jpg')
import winsound
winsound.Beep(262, 400)
winsound.Beep(262, 600)
winsound.Beep(392, 400)
winsound.Beep(349, 1000)
winsound.Beep(311, 800)
winsound.Beep(294, 800)
winsound.Beep(262, 800)
#does_optimization_help()