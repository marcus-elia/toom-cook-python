# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 08:53:58 2019

@author: Marcus

I already did all of this for standard multiplication. This is only
converting it to comvolutions.
"""
import numpy as np
import matplotlib.pyplot as plt
import time

def schoolbook(f, g, N, q):
    """ Performs convolution in R_q by doing a schoolbook multiplication
        on f and g. f and g better both have N terms."""
    print(f,g)
    product = [0]*N   # an empty list of zeros to add coefficients to
    
    for i in range(N - 1):  # -1 to do the last row separately
        for j in range(N):
            current_index = (i + j) % N
            product[current_index] += f[i]*g[j]
    
    # only mod by q at the final step
    for j in range(N):
        current_index = (- 1 + j) % N
        product[current_index] = (product[current_index] + f[-1]*g[j]) % q
    
    return product
        
# ========================
# 
# Helper Functions
#
# ========================

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


# =========================
#
# Multiplication Algorithms
#
# =========================
def karatsuba(f, g, N, q, algorithm_list):
    """ Performs a single iteration of Karatsuba multiplication and does
        schoolbook for the three smaller mults."""
    (f0, f1) = split(f,2)
    (g0, g1) = split(g,2)
    
    
    # the recursive multiplication step
    next_N = int(np.ceil(len(f) / 2))
    if not algorithm_list:
        r0 = schoolbook(f0, g0, next_N, q)
        r2 = schoolbook(f1, g1, next_N, q)
        temp_prod = schoolbook(poly_add(f0,f1), poly_add(g0,g1), next_N, q)
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r0 = karatsuba(f0, g0, next_N, q, new_list)
            r2 = karatsuba(f1, g1, next_N, q, new_list)
            temp_prod = karatsuba(poly_add(f0,f1), poly_add(g0,g1), next_N, q, new_list)
    
        elif algorithm_list[0] == '3':
            r0 = toom_3(f0, g0, next_N, q, new_list)
            r2 = toom_3(f1, g1, next_N, q, new_list)
            temp_prod = toom_3(poly_add(f0,f1), poly_add(g0,g1), next_N, q, new_list)
    
        elif algorithm_list[0] == '4':
            r0 = toom_4(f0, g0, next_N, q, new_list)
            r2 = toom_4(f1, g1, next_N, q, new_list)
            temp_prod = toom_4(poly_add(f0,f1), poly_add(g0,g1), next_N, q, new_list)
    
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    product = [0 for _ in range(N)]
    print(r0, r1, r2)
    power_of_x = 0
    for coefficient_polynomial in (r0, r1, r2):
        for i in range(len(coefficient_polynomial)):
            product[(i + power_of_x) % N] += coefficient_polynomial[i]
        power_of_x += int(np.ceil(len(f) / 2))
    return product
    
    # recombination
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    product = product[:2*len(f)-1]  
    return product

def toom_3(f, g, N, q, algorithm_list):
    """ Performs a single (optimized) iteration of Toom-3 """
     # split the polynomials into blocks
    fblocks = split(f, 3)
    gblocks = split(g, 3)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    
    next_N = int(np.ceil(len(f) / 3))
    # recursive multiplication
    if not algorithm_list:
        r = [schoolbook(f_eval[i], g_eval[i], next_N, q) for i in range(len(f_eval))]
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r = [karatsuba(f_eval[i], g_eval[i], next_N, q, new_list) for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '3':
            r = [toom_3(f_eval[i], g_eval[i], next_N, q, new_list) for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '4':
            r = [toom_4(f_eval[i], g_eval[i], N, q, new_list) for i in range(len(f_eval))]

    # this avoids matrices
    r0 = r[0]
    r4 = r[4]
    r2 = [(r[1][i] + r[2][i])//2 - r0[i] - r4[i] for i in range(len(r0))]
    r3 = [(r[3][i] - r0[i] - 4*r2[i] - 16*r4[i] + r[2][i] - r[1][i])//6 for i in range(len(r0))]
    r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] for i in range(len(r0))]
    
    #r ecombination
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
    return prod[:2*len(f)-1]
    
    
def toom_4(f, g, N, q, algorithm_list):
    """ Performs a single (optimized) iteration of Toom-4 """
    # split the polynomials into blocks
    fblocks = split(f, 4)
    gblocks = split(g, 4)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
     # recursive multiplication
    next_N = int(np.ceil(len(f) / 4))
    if not algorithm_list:
        r = {eval_list[i]:schoolbook(f_eval[i], g_eval[i], next_N, q) for i in range(len(f_eval))}
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == 'k':
            r = {eval_list[i]:karatsuba(f_eval[i], g_eval[i], next_N, q, new_list) for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '3':
            r = {eval_list[i]:toom_3(f_eval[i], g_eval[i], next_N, q, new_list) for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '4':
            r = {eval_list[i]:toom_4(f_eval[i], g_eval[i], next_N, q, new_list) for i in range(len(f_eval))}
    
    
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

def multiply(f, g, N, q, algorithm_list):
    """ This is the function we will actually plug in."""
    if not algorithm_list:
        return schoolbook(f, g, N, q)
    else:
        new_list = algorithm_list[1:]
        if algorithm_list[0] == 'k':
            product = karatsuba(f, g, N, q, new_list)
        elif algorithm_list[0] == '3':
            product = toom_3(f, g, N, q, new_list)
        elif algorithm_list[0] == '4':
            product = toom_4(f, g, N, q, new_list)
    answer = [0]*N
    for i in range(len(product)):
        answer[i % N] = (answer[i % N] + product[i]) % q
    return answer
        
f = [6,5,4,3,2,1]
g = [1,2,3,4,5,6]
print(multiply(f, g, 6, 10, ['k']))
print(schoolbook(f, g, 6, 10))