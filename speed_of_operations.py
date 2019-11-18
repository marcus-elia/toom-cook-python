# -*- coding: utf-8 -*-
"""
Created on Mon Dec 31 16:22:50 2018

@author: Marcus

This is for testing whether different operations on integers
take the same amount of time or not.
"""

import matplotlib.pyplot as plt
import time
import numpy as np

def plot_mult_vs_add(high, num_trials=10, filename=None):
    """ Makes a plot showing how long it takes on average to multiply and add
        integers up to the high value. num_trials is how many operations it
        tries at each time."""
    x_values = range(10, high, high//30)
    
    # the two y-value lists for plotting
    mult_values = []
    add_values = []
   
    for x in x_values:
        current_mult_values = []
        current_add_values = []
        for _ in range(num_trials):
            
            # two lists of random numbers for operations
            random_n1s = np.random.randint(-x, x, 10)
            random_n2s = np.random.randint(-x, x, 10)
            
            # iterate through the 100 random int pairs
            for i in range(10):            
             
                n1 = random_n1s[i]
                n2 = random_n2s[i]
                
                start_time = time.time()
                
                # add the same pair of numbers a whole bunch of times
                for _ in range(100000):
                    n1 + n2
                
                current_add_values.append(time.time() - start_time)
                
                start_time = time.time()
                
                # multiply the same pair of numbers a whole bunch of times
                for _ in range(100000):
                    n1 * n2
                
                current_mult_values.append(time.time() - start_time)
            
                    
        # put the medians in the plotting lists
        mult_values.append(np.median(current_mult_values))
        add_values.append(np.median(current_add_values))
    print(mult_values[0])
    plt.plot(x_values, add_values, lw=3, label='Addition')
    plt.plot(x_values, mult_values, lw=3, label='Multiplication')
    plt.xlabel('Integer size', size=15)
    plt.ylabel('Time (seconds)', size=15)
    plt.legend(fontsize=14)
    
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    
    plt.show()
    return

#plot_mult_vs_add(4096)
start_time = time.time()
for _ in range(100000):
    5
print(time.time() - start_time)

start_time = time.time()
for _ in range(100000):
    1000 + 1000
print(time.time() - start_time)

start_time = time.time()
for _ in range(100000):
    1000 * 1000
print(time.time() - start_time)
        

