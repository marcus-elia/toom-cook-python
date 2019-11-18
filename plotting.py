# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 22:20:31 2019

@author: Marcus
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from algorithms import *


def format_time(x):
    if x < 60:
        return str(int(x)) + " seconds"
    elif x < 3600:
        return str(int(x//60)) + "." + str(int((x-60*(x//60))/60*10)) + " minutes"
    else:
        return str(int(x//3600)) + "." + str(int((x-3600*(x//3600))/3600*10)) + " hours"



def compare_toom_algorithms(max_deg, num_trials, 
                            algorithms=(2,3,4,5,6), 
                            filename=None,
                            minsize=1):
    # the degrees (x-values)
    degrees = range(1, max_deg+1)
    
    algorithm_to_time_sum = {alg: 
        {deg:0 for deg in degrees} for alg in algorithms}
    
    # this is the outer loop so all variables will be affected equally by
    # slow outliers
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            g = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            for alg in algorithms:
                start_time = time.time()
                multiply(f, g, alg, minsize=minsize)
                total_time = time.time() - start_time
                algorithm_to_time_sum[alg][degree] += total_time
                
        # progress bar for meeeeee
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
                
    # take the average
    algorithm_to_y_values = {alg : [algorithm_to_time_sum[alg][degree]
    /num_trials for degree in degrees] for alg in algorithms}
    
    # now plot
    for alg in algorithm_to_time_sum:
        plt.plot(degrees, algorithm_to_y_values[alg], lw=3, label=str(alg))
    
    
    plt.xlabel("Degree", size=15)
    plt.ylabel("Average Time (seconds)", size=15)
    plt.legend(fontsize=13)
    plt.title("Running Time of Toom-Cook", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return

#compare_toom_algorithms(50,1000, filename='compare_toom_50.jpg')
#compare_toom_algorithms(100,1000, filename='compare_toom_100.jpg')
#compare_toom_algorithms(75, 1000, filename='compare_toom_75.jpg')
#compare_toom_algorithms(25, 1000, filename='compare_toom_25.jpg')
compare_toom_algorithms(100, 20, algorithms=([4,2],[3,3]), minsize=25)   
#f = [int(x) for x in np.random.randint(0, 2048, 50)]
#g = [int(x) for x in np.random.randint(0, 2048, 50)]
  
    
    
    
    