# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:02:24 2019

@author: Marcus


This file is for trying to use numpy's simd capabilities to implement the
vectorized index-based multiplication.
"""

import numpy as np
import time
import matplotlib.pyplot as plt

def schoolbook(f, g):
    """ Uses schoolbook multiplication to multiply f and g. Returns the
        product as a list"""
    d = len(f) + len(g) - 1
    
    # initialize a list of zeros
    product = [0]*d
    
    # distribute through all possible combinations of coefficients
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] += f[i]*g[j]
    return product

def index_based(f, g):
    t = [0]*(len(f) - 1) + f
    r = [g[0]*x for x in t]
    for i in range(1, len(f)):
        r = r[1:] + [r[0]]
        r = [r[j] + g[i]*t[j] for j in range(len(r))]
    return r

def index_based_simd(f, g):
    """ f and g are numpy arrays.
        Uses the index-based multiplication for polynomials of degree
        less than 32 using numpy for simd"""
    #t = np.array([0]*(len(f) - 1) + f)
    t = np.concatenate((np.zeros(len(f)-1, np.int16), f))
    r = g[0]*t
    for i in range(1, len(f)):
        r = np.roll(r, -1)
        r = r + g[i]*t
    return r


def format_time(x):
    """ Formats a time nicely for my progress bars when plotting"""
    if x < 60:
        return str(int(x)) + " seconds"
    elif x < 3600:
        return str(int(x//60)) + "." + str(int((x-60*(x//60))/60*10)) + " minutes"
    else:
        return str(int(x//3600)) + "." + str(int((x-3600*(x//3600))/3600*10)) + " hours"

def compare_index_based(min_deg, max_deg, num_trials, filename=None):
    
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    
    deg2schoolbook_times = {degree:0 for degree in degrees}
    deg2index_based_times = {degree:0 for degree in degrees}
    deg2simd_times = {degree:0 for degree in degrees}
    
    
    # this is the outer loop so all variables will be affected equally by
    # slow outliers
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            start_time = time.time()
            schoolbook(f, g)
            total_time = time.time() - start_time
            deg2schoolbook_times[degree] += total_time
            
            start_time = time.time()
            index_based(f, g)
            total_time = time.time() - start_time
            deg2index_based_times[degree] += total_time
            
            f = np.array(f)
            g = np.array(g)
            start_time = time.time()
            index_based_simd(f, g)
            total_time = time.time() - start_time
            deg2simd_times[degree] += total_time
            
        # progress bar for meeeeee
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
    
    # the y values to plot
    avg_schoolbook_times = []
    avg_index_based_times = []
    avg_simd_times = []
    
    for degree in degrees:
        avg_schoolbook_times.append(deg2schoolbook_times[degree] / num_trials)
        avg_index_based_times.append(deg2index_based_times[degree] / num_trials)
        avg_simd_times.append(deg2simd_times[degree] / num_trials)
        
    plt.plot(degrees, avg_schoolbook_times, label="Schoolbook", lw=3)
    plt.plot(degrees, avg_index_based_times, label="Index Based", lw=3)
    plt.plot(degrees, avg_simd_times, label="Simd", lw=3)
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Average Time (s)", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return

f = [1,2,3,4,5]
g = [6,5,4,3,2]
print(schoolbook(f,g))
print(index_based(f,g))
print(index_based_simd(f,g))            
#compare_index_based(2, 128, 300, filename='plot_for_jonathan.jpg')

f = range(1, 33)
g = [60000] + list(range(34, 65))
print([x % 2**(16) for x in schoolbook(f,g)])

