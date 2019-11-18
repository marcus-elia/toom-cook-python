# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 16:40:38 2019

@author: Marcus
Is it faster to use the temp variables in the Toom formulas?
Here is evidence:
"""


# =====================================
#
#        The stuff I always have.
#  I should put this in a separate file
#
# =====================================

import numpy as np
import time
import matplotlib.pyplot as plt

def schoolbook(f, g):
    """ Uses schoolbook multiplication to multiply f and g mod 
        m. Returns the product as a list"""
    d = len(f) + len(g) - 1
    
    # initialize a list of zeros
    product = [0]*d
    
    # distribute through all possible combinations of coefficients
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] = product[i+j] + f[i]*g[j]
    return product

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


def poly_add(f,g):
    return [f[i] + g[i] for i in range(len(f))]

def poly_subtract(f,g):
    return [f[i] - g[i] for i in range(len(f))]



# ======================================
#
#         Starting with Toom-4
# (Toom-3 might be too subtle to detect)
#
# ======================================

def plain_toom_4(f, g):
    """ Do toom4 then schoolbook"""
    n = 4    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
        
    r6 = r['infinity']
        
    r4 = [(6*r0[i] 
              - 120*r6[i] 
              + (r[2][i]+r[-2][i]) 
              - 4*(r[1][i]+r[-1][i])
              )//24 
            for i in range(len(r0))]
        
    r2 = [(r[1][i] 
              + r[-1][i])//2 
              - r0[i] 
              - r4[i] 
              - r6[i] 
              for i in range(len(r0))]
        
    r5 = [(r[3][i] 
              - 4*r[2][i] 
              + 5*r[1][i] 
              - 2*r0[i] 
              + 2*r2[i] 
              - 22*r4[i] 
              - 478*r6[i])//120 
              for i in range(len(r0))]
        
    r3 = [(r[2][i] 
              - 2*r[1][i] 
              + r0[i] 
              - 2*r2[i] 
              - 14*r4[i] 
              - 30*r5[i] 
              - 62*r6[i])//6 
             for i in range(len(r0))]
        
    r1 = [r[1][i] 
              - r0[i] 
              - r2[i] 
              - r3[i] 
              - r4[i] 
              - r5[i] 
              - r6[i] 
              for i in range(len(r0))]
   
    r_coefs = (r0, r1, r2, r3, r4, r5, r6)
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def temp_toom_4(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    
    L = range(len(r0)) # to save time on for-loops
        
    r6 = r['infinity']

    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r6[i]
          for i in L]
    
    r4 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 64*r6[i]) // 4 - T1[i]) // 3
          for i in L]
    
    r2 = [T1[i] - r4[i] for i in L]
    
    T2 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    
    T3 = [((r[2][i] - r[-2][i]) // 4 - T2[i]) // 3 for i in L]
    
    T4 = [((r[3][i] - r0[i] - 9*r2[i] - 81*r4[i] - 729*r6[i]) // 3
           - T2[i]) // 8
         for i in L]

    r5 = [(T4[i] - T3[i]) // 5 for i in L]
    
    r3 = [T3[i] - 5*r5[i] for i in L]
    
    r1 = [T2[i] - r3[i] - r5[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6)

    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def plain_toom_3(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    r0 = r[0]
        
    r4 = r['infinity']
        
    r2 = [(r[1][i] + r[-1][i])//2 
          - r0[i] 
          - r4[i] 
          for i in range(len(r0))]

    r3 = [(r[2][i] 
          - r0[i] 
          - 4*r2[i] 
          - 16*r4[i] 
          + r[-1][i] 
          - r[1][i])//6 
         for i in range(len(r0))]
        
    r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] for i in range(len(r0))]
    
    r_coefs = (r0, r1, r2, r3, r4)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def temp_toom_3(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    
    L = range(len(r0)) # to save time on for-loops
        
    r4 = r['infinity']
     
    r2 = [(r[1][i] + r[-1][i])//2 
              - r0[i] 
              - r4[i] 
              for i in L]
    
    T = [(r[1][i] - r[-1][i]) // 2 for i in L]
    
    r3 = [((r[2][i] - r0[i] - 4*r2[i] - 16*r4[i]) // 2 - T[i]) // 3
           for i in L]
    
    r1 = [T[i] - r3[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def plot_times_toom_n(n, min_deg, max_deg, num_trials, filename=None):
    degrees = range(min_deg, max_deg+1)
    plain_time = {degree:0 for degree in degrees}
    temp_time = {degree:0 for degree in degrees}
    
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            if n == 3:
                start_time = time.time()
                plain_toom_3(f, g)
                plain_time[degree] += time.time() - start_time
            
                start_time = time.time()
                temp_toom_3(f, g)
                temp_time[degree] += time.time() - start_time
            if n == 4:
                start_time = time.time()
                plain_toom_4(f, g)
                plain_time[degree] += time.time() - start_time
            
                start_time = time.time()
                temp_toom_4(f, g)
                temp_time[degree] += time.time() - start_time
            
        # progress bar
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
    
    plain_time = [plain_time[degree]/num_trials for degree in degrees]
    temp_time = [temp_time[degree]/num_trials for degree in degrees]
    
    plt.plot(degrees, plain_time, label="Plain Toom-{}".format(n))
    plt.plot(degrees, temp_time, label="Temp Toom-{}".format(n))
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Time (seconds)", size=15)
    plt.title("Average Time for Polynomial Multiplication", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return
    
    
#plot_times_toom_4(80, 120, 2000, filename='is_temp_4_faster.jpg')    
plot_times_toom_n(3, 61, 101, 10000, filename='is_temp_3_faster.jpg')    
    
    
    
    
    
    
    
    
    
    
    