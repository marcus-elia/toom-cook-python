# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 08:37:59 2019

@author: Marcus
"""

from algorithms import *
import numpy as np
import matplotlib.pyplot as plt
import time

# -----------------------------
# =============================
#
#      Helper Functions
#
# =============================
# -----------------------------

def format_time(x):
    if x < 60:
        return str(int(x)) + " seconds"
    elif x < 3600:
        return str(int(x//60)) + "." + str(int((x-60*(x//60))/60*10)) + " minutes"
    else:
        return str(int(x//3600)) + "." + str(int((x-3600*(x//3600))/3600*10)) + " hours"
    
def list_to_string(algorithm_list):
    """ Converts a list of ints into a string with a dash between
        each character"""
    if not algorithm_list:
        return "[]"
    if len(algorithm_list) == 1:
        return str(algorithm_list[0])
    s = str(algorithm_list[0])
    for n in algorithm_list[1:-1]:
        s += '-'
        s += str(n)
    return s + '-' + str(algorithm_list[-1])

# -----------------------------
# =============================
#
#         Plotting
#
# =============================
# -----------------------------

def compare_algorithms(min_deg, max_deg, num_trials, 
                            algorithms, 
                            filename=None):
    algorithms = [tuple(algorithm_list) for algorithm_list in algorithms]
    
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    
    # a dictionary mapping each algorithm list to a dictionary mapping each
    # degree to the time sum for it
    algorithm_to_time_sum = {algorithm_list: {degree:0 for degree in degrees} 
                             for algorithm_list in algorithms}
    
    # this is the outer loop so all variables will be affected equally by
    # slow outliers
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            for algorithm_list in algorithms:
                start_time = time.time()
                multiply(f, g, algorithm_list)
                total_time = time.time() - start_time
                algorithm_to_time_sum[algorithm_list][degree] += total_time
                
        # progress bar for meeeeee
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
                
    # take the average
    algorithm_list_to_y_values = {algorithm_list :
        [algorithm_to_time_sum[algorithm_list][degree]/num_trials for 
         degree in degrees] for algorithm_list in algorithms}
    
    # now plot
    for algorithm_list in algorithm_to_time_sum:
        plt.plot(degrees, algorithm_list_to_y_values[algorithm_list], 
                 lw=3, label=list_to_string(algorithm_list))

    plt.xlabel("Degree", size=15)
    plt.ylabel("Average Time (seconds)", size=15)
    plt.legend(fontsize=13)
    plt.title("Running Time of Multiplication Algorithms", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return

# -----------------------------
# =============================
#
#          Timing
#
# =============================
# -----------------------------

def print_algorithm_times(degree, num_trials, algorithms):
    """ Makes a dictionary mapping each algorithm to the average amount
        of time it takes"""
        
    # dictionaries don't like lists
    algorithms = [tuple(algorithm) for algorithm in algorithms]
    algorithm_to_time = {algorithm:0 for algorithm in algorithms}

    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        
        # random polynomials
        f = [int(x) for x in np.random.randint(0, 2048, degree)]
        g = [int(x) for x in np.random.randint(0, 2048, degree)]
        
        # measure the time for each algorithms
        for algorithm in algorithms:
            start_time = time.time()
            multiply(f, g, algorithm)
            total_time = time.time() - start_time
            algorithm_to_time[algorithm] += total_time

        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
    
    

    for algorithm in algorithm_to_time:
        print("{}: {}".format(algorithm,
              algorithm_to_time[algorithm]/num_trials))
    return

compare_algorithms(735, 750, 200, [[2,2,2,2,2],[4,3,2]], filename='2vs432.jpg')
#compare_algorithms(50, 80, 500, [[3],[2,2]])
#print_algorithm_times(742, 1000, [[6,4],[8,3], [5,5], [12,2]])
