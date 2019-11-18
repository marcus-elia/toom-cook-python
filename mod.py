# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:24:02 2019

@author: Marcus
This is the same as compound multiplication algorithms, except mod. And the
actual testing is done with one trinary polynomial multiplied by one random
polynomial.
"""

import numpy as np
import matplotlib.pyplot as plt
import time


# ==========================================
#
# Helper functions for all of the algorithms
#
# ==========================================

def split(f, num_blocks):
    """ Splits the list f into num_blocks different blocks of equal size
        If it doesn't divide evenly, we put zeros on the end of the last
        block."""
    blocks = []
    copy_f = list(f)  # copy f so we don't ruin it!!!!!!!!
    while len(copy_f) % num_blocks != 0:
        copy_f.append(0)
    block_length = len(copy_f) // num_blocks
    index = 0
    while index + block_length < len(copy_f):
        blocks.append(copy_f[index:index+block_length])
        index += block_length
    blocks.append(copy_f[index:])
    return blocks

def evaluate_blocks(blocks, value):
    """ blocks is a list of lists, each list is the coefficients of a
        polynomial. But each list a coefficient. For example, if blocks is
        [[1,2],[3,4],[5,6]] and value is -2, we return
        [1,2] + [-6,-8] + [20,24] = [15, 18].  If the value is infinity,
        we return the leading coefficient."""
        
    if value == 'infinity':
        return blocks[-1]
    
    # initialize an empty list of the right length
    answer = [0]*len(blocks[0])
    
    coefficient = 1
    for i in range(len(blocks)):
        for j in range(len(blocks[0])):
            answer[j] += coefficient*blocks[i][j]
        coefficient *= value    # multiply to make powers of value
    return answer

def evaluate_blocks_list(blocks, values):
    """ Evaluates the blocks on a list of values, and returns a list"""
    answer = []
    for value in values:
        answer.append(evaluate_blocks(blocks, value))
    return answer

def poly_add(f,g):
    return [f[i] + g[i] for i in range(len(f))]

def poly_subtract(f,g):
    return [f[i] - g[i] for i in range(len(f))]

def list_to_string(algorithm_list):
    """ Converts a list of characters into a string with a dash between
        each character"""
    if not algorithm_list:
        return "[]"
    if len(algorithm_list) == 1:
        return algorithm_list[0]
    s = algorithm_list[0]
    for char in algorithm_list[1:-1]:
        s += '-'
        s += char
    return s + '-' + algorithm_list[-1]

def format_time(x):
    if x < 60:
        return str(int(x)) + " seconds"
    elif x < 3600:
        return str(int(x//60)) + "." + str(int((x-60*(x//60))/60*10)) + " minutes"
    else:
        return str(int(x//3600)) + "." + str(int((x-3600*(x//3600))/3600*10)) + " hours"

# Multiplication Algorithms

def schoolbook_multiply(f, g, q):
    """ Just plain multiplies f and g by distributing. Returns the
        the answer mod q."""
    d = len(f) + len(g) - 1
    
    # create a list of zeros
    product = [0]*d
    
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] =  (product[i + j] + f[i]*g[j]) % q
            
    return product

# In all of these algorithms, we input a list of algorithms to use next.
# For example ['k', '4', '3'] means do Karatsuba next, then Toom 4, then
# Toom 3, and then just schoolbook.

def karatsuba(f,g, algorithm_list, q):
    """ Performs a single iteration of Karatsuba multiplication and does
        schoolbook for the three smaller mults."""
    (f0, f1) = split(f,2)
    (g0, g1) = split(g,2)
    
    
    # the recursive multiplication step
    if not algorithm_list:
        r0 = schoolbook_multiply(f0, g0, q)
        r2 = schoolbook_multiply(f1, g1, q)
        temp_prod = schoolbook_multiply(poly_add(f0,f1), poly_add(g0,g1), q)
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r0 = karatsuba(f0, g0, new_list, q)
            r2 = karatsuba(f1, g1, new_list, q)
            temp_prod = karatsuba(poly_add(f0,f1), 
                                  poly_add(g0,g1), 
                                  new_list, 
                                  q)
    
        elif algorithm_list[0] == '3':
            r0 = toom_3(f0, g0, new_list, q)
            r2 = toom_3(f1, g1, new_list, q)
            temp_prod = toom_3(poly_add(f0,f1), poly_add(g0,g1), new_list, q)
    
        elif algorithm_list[0] == '4':
            r0 = toom_4(f0, g0, new_list, q)
            r2 = toom_4(f1, g1, new_list, q)
            temp_prod = toom_4(poly_add(f0,f1), poly_add(g0,g1), new_list, q)
    
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    
    # this is the most confusing thing ever
    if len(f) % 2 == 0:
        k = len(f)//2
        part1 = r0[:k] 
        part2 = [r0[k+i] + r1[i] for i in range(k-1)] 
        part3 = [r1[k-1]]
        part4 = [r1[k+i] + r2[i] for i in range(k-1)]
        part5 = r2[k-1:]
        return part1 + part2 + part3 + part4 + part5
    else:
        k = len(f)//2
        part1 = r0[:k+1]
        part2 = [r0[k+1+i] + r1[i] for i in range(k)]
        part3 = [r1[k]]
        part4 = [r1[k+i+1] + r2[i] for i in range(0,k)]
        part5 = r2[k:4*k+1 - 2*k-2]
        return part1 + part2 + part3 + part4 + part5   

def toom_3(f,g, algorithm_list, q):
    """ Performs a single (optimized) iteration of Toom-3 """
     # split the polynomials into blocks
    fblocks = split(f, 3)
    gblocks = split(g, 3)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    
    # recursive multiplication
    if not algorithm_list:
        r = [schoolbook_multiply(f_eval[i], g_eval[i], q) for i in range(len(f_eval))]
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r = [karatsuba(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '3':
            r = [toom_3(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '4':
            r = [toom_4(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))]
    print(r)
    # this avoids matrices
    r0 = r[0]
    r4 = r[4]
    r2 = [(r[1][i] + r[2][i])//2 - r0[i] - r4[i] for i in range(len(r0))]
    r3 = [(r[3][i] - r0[i] - 4*r2[i] - 16*r4[i] + r[2][i] - r[1][i])//6 for i in range(len(r0))]
    r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] for i in range(len(r0))]
    print(r0, r1, r2, r3, r4)
    
    # recombination
    k = int(np.ceil(len(f) / 3))
    prod = r0[:k]
    prod = prod + [r0[k+i] + r1[i] for i in range(k-1)]
    prod = prod + [r1[k-1]]
    prod = prod + [r1[k+i] + r2[i] for i in range(k-1)]
    prod = prod + [r2[k-1]]
    prod = prod + [r2[k+i] + r3[i] for i in range(k-1)]
    prod = prod + [r3[k-1]]
    prod = prod + [r3[k+i] + r4[i] for i in range(k-1)]
    prod = prod + r4[k-1:]
    
    # I couldn't quite figure out how to end the line above this one (it was
    # only working in 0 mod 3 and 2 mod 3), so I'm just chopping off the 
    # excess 0's here, if there are any.
    return prod[:2*len(f) - 1]

def toom_4(f ,g, algorithm_list, q):
    """ Performs a single (optimized) iteration of Toom-4 """
    # split the polynomials into blocks
    fblocks = split(f, 4)
    gblocks = split(g, 4)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
     # recursive multiplication
    if not algorithm_list:
        r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i], q) for i in range(len(f_eval))}
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r = {eval_list[i]:karatsuba(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '3':
            r = {eval_list[i]:toom_3(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '4':
            r = {eval_list[i]:toom_4(f_eval[i], g_eval[i], new_list, q) for i in range(len(f_eval))}
    
    
    r0 = r[0]
    r6 = r['infinity']
    r4 = [(6*r0[i] - 120*r6[i] + (r[2][i]+r[-2][i]) - 4*(r[1][i]+r[-1][i]))//24 for i in range(len(r0))]
    r2 = [(r[1][i] + r[-1][i])//2 - r0[i] - r4[i] - r6[i] for i in range(len(r0))]
    r5 = [(r[3][i] - 4*r[2][i] + 5*r[1][i] - 2*r0[i] + 2*r2[i] - 22*r4[i] - 478*r6[i])//120 for i in range(len(r0))]
    r3 = [(r[2][i] - 2*r[1][i] + r0[i] -2*r2[i] - 14*r4[i] -30*r5[i] - 62*r6[i])//6 for i in range(len(r0))]
    r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] - r5[i] - r6[i] for i in range(len(r0))]
    
    # need these in a tuple to automate recombination
    r_coefs = (r0, r1, r2, r3, r4, r5, r6)
    
    # recombination
    k = int(np.ceil(len(f) / 4))
    prod = r0[:k]
    for j in range(1,6):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]
    prod = prod + [r_coefs[5][k+i] + r_coefs[6][i] for i in range(k-1)]
    prod = prod + r_coefs[6][k-1:]
    
    return prod[:2*len(f)-1]

def multiply(f, g, algorithm_list, q):
    """ This is the function we will actually plug in."""
    if not algorithm_list:
        return schoolbook_multiply(f, g, q)
    else:
        new_list = algorithm_list[1:]
        if algorithm_list[0] == 'k':
            product = karatsuba(f, g, new_list, q)
        elif algorithm_list[0] == '3':
            product = toom_3(f, g, new_list, q)
        elif algorithm_list[0] == '4':
            product = toom_4(f, g, new_list, q)
    
    return [x % q for x in product]

def plot_multiple_algorithms(min_deg, max_deg, num_trials, algorithms, filename=None):
    """ Specify the minimum and maximum degree for our plot, and the number of
        trials at each degree to average over. Then specify a list of
        algorithms, where each algorithm is like ['k', '3', '4']."""
    algorithms = [tuple(algorithm_list) for algorithm_list in algorithms]
    
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    
    # a dictionary mapping each algorithm list to a dictionary mapping each
    # degree to the time sum for it
    algorithm_to_time_sum = {algorithm_list: {degree:0 for degree in degrees} for algorithm_list in algorithms}
    
    # this is the outer loop so all variables will be affected equally by
    # slow outliers
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(-2048, 2048, degree)]
            g = [int(x) for x in np.random.randint(-2048, 2048, degree)]
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
    algorithm_list_to_y_values = {algorithm_list : [algorithm_to_time_sum[algorithm_list][degree]/num_trials for degree in degrees] for algorithm_list in algorithms}
    
    # now plot
    for algorithm_list in algorithm_to_time_sum:
        plt.plot(degrees, algorithm_list_to_y_values[algorithm_list], lw=3, label=list_to_string(algorithm_list))

    plt.xlabel("Degree", size=15)
    plt.ylabel("Average Time (seconds)", size=15)
    plt.legend(fontsize=13)
    plt.title("Running Time of Multiplication Algorithms", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return

f = [1,2,3,4,5,6]
g = [6,5,4,3,2,1]
print(schoolbook_multiply(f,g, 10))
print(multiply(f, g, ['3'], 10))