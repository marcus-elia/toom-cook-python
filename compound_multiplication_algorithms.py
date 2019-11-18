# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 09:22:55 2019

@author: Marcus

This file is for implementing different combinations of Toom-Cook
algorithms for the purpose of multiplying large (742 degree) polynomials.
"""
import numpy as np
import matplotlib.pyplot as plt
import time


def num_int_operations(degree, algorithm_list):
    """ degree is the degree of the polynomials we want to multiply.
        algorithm_list is a list of characters from the set 
        {2, 3, 4} to specify Karatsuba, Toom-3, and Toom-4. After 
        obeying the list, schoolbook is used. For example,
        num_int_operations(768, ['2', '4', '3']) = 327452."""
    
    # the base case: schoolbook multiplication. This takes 2n^2 + 2n + 1
    # integer operations
    if not algorithm_list:
        return 5*(degree**2) + 12*degree + 8
    
    if algorithm_list[0] == '2':
        extras = 46*degree + 45
        return 3*num_int_operations(degree//2, algorithm_list[1:]) + extras
    
    if algorithm_list[0] == '3':
        extras = 272*degree/3 + 352/3
        return 5*num_int_operations(degree//3, algorithm_list[1:]) + extras
    
    if algorithm_list[0] == '4':
        extras = 269*degree/2 + 523/2
        return 7*num_int_operations(degree//4, algorithm_list[1:]) + extras
    
#print(num_int_operations(768, ['k', '3', '4']))
#print(num_int_operations(742,['4','4','4']))

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

def schoolbook_multiply(f,g):
    """ Just plain multiplies f and g by distributing"""
    d = len(f) + len(g) - 1
    
    # create a list of zeros
    product = [0]*d
    
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] += f[i]*g[j]
            
    return product

# In all of these algorithms, we input a list of algorithms to use next.
# For example ['k', '4', '3'] means do Karatsuba next, then Toom 4, then
# Toom 3, and then just schoolbook.

def karatsuba(f,g, algorithm_list):
    """ Performs a single iteration of Karatsuba multiplication and does
        schoolbook for the three smaller mults."""
    (f0, f1) = split(f,2)
    (g0, g1) = split(g,2)
    
    
    # the recursive multiplication step
    if not algorithm_list:
        r0 = schoolbook_multiply(f0, g0)
        r2 = schoolbook_multiply(f1, g1)
        temp_prod = schoolbook_multiply(poly_add(f0,f1), poly_add(g0,g1))
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == '2':
            r0 = karatsuba(f0, g0, new_list)
            r2 = karatsuba(f1, g1, new_list)
            temp_prod = karatsuba(poly_add(f0,f1), poly_add(g0,g1), new_list)
    
        elif algorithm_list[0] == '3':
            r0 = toom_3(f0, g0, new_list)
            r2 = toom_3(f1, g1, new_list)
            temp_prod = toom_3(poly_add(f0,f1), poly_add(g0,g1), new_list)
    
        elif algorithm_list[0] == '4':
            r0 = toom_4(f0, g0, new_list)
            r2 = toom_4(f1, g1, new_list)
            temp_prod = toom_4(poly_add(f0,f1), poly_add(g0,g1), new_list)
        
        elif algorithm_list[0] == '5':
            r0 = toom_5(f0, g0, new_list)
            r2 = toom_5(f1, g1, new_list)
            temp_prod = toom_5(poly_add(f0,f1), poly_add(g0,g1), new_list)
            
        elif algorithm_list[0] == '6':
            r0 = toom_6(f0, g0, new_list)
            r2 = toom_6(f1, g1, new_list)
            temp_prod = toom_6(poly_add(f0,f1), poly_add(g0,g1), new_list)
    
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

def toom_3(f,g, algorithm_list):
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
        r = [schoolbook_multiply(f_eval[i], g_eval[i])
        for i in range(len(f_eval))]
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == '2':
            r = [karatsuba(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '3':
            r = [toom_3(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '4':
            r = [toom_4(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '5':
            r = [toom_5(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))]
            
        elif algorithm_list[0] == '6':
            r = [toom_6(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))]

    # this avoids matrices
    r0 = r[0]
    r4 = r[4]
    r2 = [(r[1][i] + r[2][i])//2 - r0[i] - r4[i] for i in range(len(r0))]
    r3 = [(r[3][i] - r0[i] - 4*r2[i] - 16*r4[i] + r[2][i] - r[1][i])//6 for i in range(len(r0))]
    r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] for i in range(len(r0))]
    
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

def toom_4(f ,g, algorithm_list):
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
        r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i])
        for i in range(len(f_eval))}
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == '2':
            r = {eval_list[i]:karatsuba(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '3':
            r = {eval_list[i]:toom_3(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '4':
            r = {eval_list[i]:toom_4(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '5':
            r = {eval_list[i]:toom_5(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '6':
            r = {eval_list[i]:toom_6(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
    
    
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


def toom_5(f ,g, algorithm_list):
    """ Performs a single (optimized) iteration of Toom-5 """
    # split the polynomials into blocks
    fblocks = split(f, 5)
    gblocks = split(g, 5)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
     # recursive multiplication
    if not algorithm_list:
        r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i])
        for i in range(len(f_eval))}
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == '2':
            r = {eval_list[i]:karatsuba(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '3':
            r = {eval_list[i]:toom_3(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '4':
            r = {eval_list[i]:toom_4(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '5':
            r = {eval_list[i]:toom_5(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '6':
            r = {eval_list[i]:toom_6(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
    
    
    r0 = r[0]
    
    r8 = r['infinity']
    
    r6 = [(r[3][i] + r[-3][i]
           - 6*(r[2][i] + r[-2][i]) 
           + 15*(r[1][i] + r[-1][i]) 
           - 20*r0[i] 
           - 10080*r8[i])//720 
           for i in range(len(r0))]
    
    r4 = [((r[2][i] + r[-2][i])//2 
           - 2*(r[1][i] + r[-1][i]) 
          + 3*r0[i] 
          - 60*r6[i] 
          - 252*r8[i])//12 
          for i in range(len(r0))]
    
    r2 = [(r[1][i] + r[-1][i])//2 
          - r0[i] 
          - r4[i] 
          - r6[i] 
          - r8[i] 
          for i in range(len(r0))]
    
    r7 = [(r[4][i] 
           - 14*r[1][i] 
           + 14*r[2][i] 
           - 6*r[3][i] 
           + 5*r0[i] 
           - 4*r2[i] 
           + 20*r4[i] 
           - 604*r6[i] 
           - 29740*r8[i]) // 5040 
           for i in range(len(r0))]
    
    r5 = [(r[3][i] 
           + 5*r[1][i] 
           - 4*r[2][i] 
           - 2*r0[i] 
           + 2*r2[i] 
           - 22*r4[i] 
           - 478*r6[i] 
           - 1680*r7[i] 
           - 5542*r8[i]) // 120 
           for i in range(len(r0))]
    
    r3 = [(r[2][i] 
           - 2*r[1][i] 
           + r0[i] 
           - 2*r2[i] 
           - 14*r4[i] 
           - 30*r5[i] 
           - 62*r6[i] 
           -126*r7[i] 
           - 254*r8[i]) // 6 
           for i in range(len(r0))]
    
    r1 = [r[1][i]
           - r0[i] 
           - r2[i] 
           - r3[i] 
           - r4[i] 
           - r5[i] 
           - r6[i] 
           - r7[i] 
           - r8[i] 
           for i in range(len(r0))]
    
    
    # need these in a tuple to automate recombination
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
       
    # recombination
    k = int(np.ceil(len(f) / 5))
    prod = r0[:k]
    for j in range(1,8):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]
    prod = prod + [r_coefs[7][k+i] + r_coefs[8][i] for i in range(k-1)]
    prod = prod + r_coefs[8][k-1:]
    
    return prod[:2*len(f)-1]


def toom_6(f ,g, algorithm_list):
    """ Performs a single (optimized) iteration of Toom-6 """
    # split the polynomials into blocks
    fblocks = split(f, 6)
    gblocks = split(g, 6)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
     # recursive multiplication
    if not algorithm_list:
        r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i])
        for i in range(len(f_eval))}
    else:
        new_list = algorithm_list[1:]
        
        if algorithm_list[0] == '2':
            r = {eval_list[i]:karatsuba(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '3':
            r = {eval_list[i]:toom_3(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '4':
            r = {eval_list[i]:toom_4(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
        elif algorithm_list[0] == '5':
            r = {eval_list[i]:toom_5(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
        
        elif algorithm_list[0] == '6':
            r = {eval_list[i]:toom_6(f_eval[i], g_eval[i], new_list)
            for i in range(len(f_eval))}
            
    r0 = r[0]
    
    r10 = r['infinity']
    
    r8 = [((r[4][i] + r[-4][i])
           - 8*(r[3][i] + r[-3][i]) 
           + 28*(r[2][i] + r[-2][i]) 
           - 56*(r[1][i] + r[-1][i]) 
           + 70*r0[i] 
           - 1209600*r10[i])//40320 
        for i in range(len(r0))]
    
    r6 = [((r[3][i] + r[-3][i]) 
          - 6*(r[2][i] + r[-2][i]) 
          + 15*(r[1][i] + r[-1][i]) 
          - 20*r0[i] 
          - 10080*r8[i] 
          - 105840*r10[i])//720 
        for i in range(len(r0))]
    
    r4 = [((r[2][i] + r[-2][i]) 
          - 4*(r[1][i] + r[-1][i]) 
          + 6*r0[i] 
          - 120*r6[i] 
          - 504*r8[i] 
          - 2040*r10[i])//24 
        for i in range(len(r0))]
    
    r2 = [(r[1][i] + r[-1][i])//2 
          - r0[i] 
          - r4[i] 
          - r6[i] 
          - r8[i] 
          - r10[i] 
          for i in range(len(r0))]
    
    r9 = [(5*r[5][i]
           - 40*r[4][i] 
           + 135*r[3][i] 
           - 240*r[2][i] 
           + 210*r[1][i]
           - 70*r0[i] 
           + 50*r2[i] 
           - 190*r4[i] 
           + 2450*r6[i] 
           - 156190*r8[i] 
           - 14611150*r10[i])//1814400 
         for i in range(len(r0))]
    
    r7 = [(r[5][i] 
           - 2*r[4][i] 
           - 9*r[3][i] 
           + 36*r[2][i] 
           - 42*r[1][i] 
           + 16*r0[i] 
           - 14*r2[i] 
           + 82*r4[i] 
           - 3134*r6[i] 
           - 209678*r8[i] 
           - 1270080*r9[i] 
           - 7173854*r10[i])//30240 
          for i in range(len(r0))]
    
    r5 = [(r[3][i] 
           - 4*r[2][i] 
           + 5*r[1][i] 
           - 2*r0[i] 
           + 2*r2[i] 
           - 22*r4[i] 
           - 478*r6[i] 
           - 1680*r7[i] 
           - 5542*r8[i] 
           - 17640*r9[i] 
           - 54958*r10[i])//120 
         for i in range(len(r0))]
    
    r3 = [(r[2][i] 
           - 2*r[1][i] 
           + r0[i] 
           - 2*r2[i] 
           - 14*r4[i] 
           - 30*r5[i] 
           - 62*r6[i] 
           - 126*r7[i] 
           - 254*r8[i] 
           - 510*r9[i] 
           - 1022*r10[i])//6 
         for i in range(len(r0))]
    
    r1 = [r[1][i] 
        - r0[i] 
        - r2[i] 
        - r3[i] 
        - r4[i] 
        - r5[i] 
        - r6[i] 
        - r7[i] 
        - r8[i] 
        - r9[i] 
        - r10[i] 
        for i in range(len(r0))]
    
    # need these in a tuple to automate recombination
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    
    # recombination
    k = int(np.ceil(len(f) / 6))
    prod = r0[:k]
    for j in range(1,10):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]
    prod = prod + [r_coefs[9][k+i] + r_coefs[10][i] for i in range(k-1)]
    prod = prod + r_coefs[10][k-1:]
    
    return prod[:2*len(f)-1]


# ================================
#
# This is where it's all put together
#
# =================================

def multiply(f, g, algorithm_list, q=None):
    """ This is the function we will actually plug in. If q is specified,
        then we mod the coefficients by q at the very last step."""
    if not algorithm_list:
        return schoolbook_multiply(f,g)
    else:
        new_list = algorithm_list[1:]
        if algorithm_list[0] == '2':
            product = karatsuba(f, g, new_list)
        elif algorithm_list[0] == '3':
            product = toom_3(f, g, new_list)
        elif algorithm_list[0] == '4':
            product = toom_4(f, g, new_list)
        elif algorithm_list[0] == '5':
            product = toom_5(f, g, new_list)
        elif algorithm_list[0] == '6':
            product = toom_6(f, g, new_list)
    if not q:
        return product
    else:
        return [x % q for x in product]
    
    
# =====================================    

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

def plot_multiple_algorithms_NTRU(min_deg, max_deg, num_trials, 
                                    algorithms, filename=None):
    """ Specify the minimum and maximum degree for our plot, and the number of
        trials at each degree to average over. Then specify a list of
        algorithms, where each algorithm is like ['k', '3', '4']. This
        happens mod 2048, and we multiply trinomials by random polynomials."""
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
            
            # make a random trinomial
            trinom1 = [1]*(degree//3) 
            trinom2 = [-1]*(degree//3) 
            trinom3 = [0]*(degree - degree//3 - degree//3)
            f = np.random.permutation(trinom1 + trinom2 + trinom3)
            
            # make a random polynomial mod 2048
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            for algorithm_list in algorithms:
                start_time = time.time()
                multiply(f, g, algorithm_list, q=2048)
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
    plt.title("Running Time of Multiplication Algorithms for NTRU", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return


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


#print_algorithm_times(742, 2000, [['6','4'], ['5','2','2'], ['6','5'],
#                                 ['6','2','2'], ['3','3','3'], ['4','4','3'],
#                                 ['3','4','3'], ['4','4','4'], [],
#                                 ['2','2','2','2','2'], ['5','3','2'],
#                                 ['5','5'], ['4','3','2'], ['2','4','3'],
#                                 ['4','4','2'], ['3','4','2'], ['5','4'],
#                                 ['5','4','2'], ['6','3','2'], ['2','3','4'],
#                                 ['4', '3', '3'], ['6','6']])

#print_algorithm_times(442, 2000, [['5','3'], ['4','4'], ['4','2','2'],
#                                 ['6','3'], ['3','3','2'], ['5','4']])

#plot_multiple_algorithms(735, 750, 200, [['6', '5'], ['4','4', '2']])

if False:
    import winsound
    winsound.Beep(262, 400)
    winsound.Beep(262, 600)
    winsound.Beep(392, 400)
    winsound.Beep(349, 1000)
    winsound.Beep(311, 800)
    winsound.Beep(294, 800)
    winsound.Beep(262, 800)