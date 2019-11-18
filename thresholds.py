# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 07:24:43 2019

@author: Marcus

This file is for very precisely computing the thresholds of various algorithms
"""
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

# =======================================
#
#          Single Algorithms
#
# =======================================

def toom2(f, g):
    (f0, f1) = split(f, 2)
    (g0, g1) = split(g, 2)
    
    r0 = schoolbook(f0, g0)
    r2 = schoolbook(f1, g1)
    temp_prod = schoolbook(poly_add(f0, f1), poly_add(g0, g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1] 

def toom3(f, g):
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

def toom4(f, g):
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


def toom5(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    
    r8 = r['infinity']
    
    L = range(len(r0)) # to save time on for-loops
    
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    
    r4 = [T2[i] - 5*r6[i] for i in L]
    
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    
    r5 = [T8[i] - 14*r7[i] for i in L]
    
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]
 
    
def toom6(f, g):
    n = 6
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r10 = r['infinity']
    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i]
          for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 1024*r10[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 59049*r10[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 1048576*r10[i])//16 - T1[i])//15
          for i in L]
    T5 = [(T3[i] - T2[i])//5 for i in L]
    T6 = [(T4[i] - T3[i])//7 for i in L]
    r8 = [(T6[i] - T5[i])//12 for i in L]
    r6 = [T5[i] - 14*r8[i] for i in L]
    r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in L]
    
    T7 = [(r[1][i] - r[-1][i])//2 for i in L]
    T8 = [((r[2][i] - r[-2][i])//4 - T7[i])//3 for i in L]
    T9 = [((r[3][i] - r[-3][i])//6 - T7[i])//8 for i in L]
    T10 = [((r[4][i] - r[-4][i])//8 - T7[i])//15 for i in L]
    T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] - 15625*r6[i] 
    - 390625*r8[i] - 9765625*r10[i])//5 - T7[i])//24 for i in L]
    T12 = [(T9[i] - T8[i])//5 for i in L]
    T13 = [(T10[i] - T9[i])//7 for i in L]
    T14 = [(T11[i] - T10[i])//9 for i in L]
    T15 = [(T13[i] - T12[i])//12 for i in L]
    T16 = [(T14[i] - T13[i])//16 for i in L]
    r9 = [(T16[i] - T15[i])//21 for i in L]
    r7 = [T15[i] - 30*r9[i] for i in L]
    r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in L]
    r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in L]
    r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom7(f, g):
    n = 7
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
    L = range(len(r0))
    r12 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r12[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 4096*r12[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 531441*r12[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 16777216*r12[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 244140625*r12[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(T3[i] - T2[i]) // 5 for i in L]
    T7 = [(T4[i] - T3[i]) // 7 for i in L]
    T8 = [(T5[i] - T4[i]) // 9 for i in L]
    T9 = [(T7[i] - T6[i]) // 12 for i in L]
    T10 = [(T8[i] - T7[i]) // 16 for i in L]
    r10 = [(T10[i] - T9[i]) // 21 for i in L]
    r8 = [T9[i] - 30*r10[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] for i in L]
    T11 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T12 = [((r[2][i] - r[-2][i]) // 4 - T11[i]) // 3 for i in L]
    T13 = [((r[3][i] - r[-3][i])//6 - T11[i])//8 for i in L]
    T14 = [((r[4][i] - r[-4][i])//8 - T11[i])//15 for i in L]
    T15 = [((r[5][i] - r[-5][i])//10 - T11[i])//24 for i in L]
    T16 = [((r[6][i] - r0[i] - 36*r2[i] - 1296*r4[i] - 46656*r6[i] 
           - 1679616*r8[i] - 60466176*r10[i] - 2176782336*r12[i])//6 
           - T11[i])//35 for i in L]
    T17 = [(T13[i] - T12[i]) // 5 for i in L]
    T18 = [(T14[i] - T13[i]) // 7 for i in L]
    T19 = [(T15[i] - T14[i]) // 9 for i in L]
    T20 = [(T16[i] - T15[i]) // 11 for i in L]
    T21 = [(T18[i] - T17[i]) // 12 for i in L]
    T22 = [(T19[i] - T18[i]) // 16 for i in L]
    T23 = [(T20[i] - T19[i]) // 20 for i in L]
    T24 = [(T22[i] - T21[i]) // 21 for i in L]
    T25 = [(T23[i] - T22[i]) // 27 for i in L]
    r11 = [(T25[i] - T24[i]) // 32 for i in L]
    r9 = [T24[i] - 55*r11[i] for i in L]
    r7 = [T21[i] - 30*r9[i] - 627*r11[i] for i in L]
    r5 = [T17[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] for i in L]
    r3 = [T12[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i]
          for i in L]
    r1 = [T11[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom8(f, g):
    n = 8
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, -6, 7, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r14 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r14[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 16384*r14[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 4782969*r14[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 268435456*r14[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 6103515625*r14[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(((r[6][i] + r[-6][i])//2 - r0[i] - 78364164096*r14[i])//36 
           - T1[i])//35 for i in L]
    T7 = [(T3[i] - T2[i]) // 5 for i in L]
    T8 = [(T4[i] - T3[i]) // 7 for i in L]
    T9 = [(T5[i] - T4[i]) // 9 for i in L]
    T10 = [(T6[i] - T5[i]) // 11 for i in L]
    T11 = [(T8[i] - T7[i]) // 12 for i in L]
    T12 = [(T9[i] - T8[i]) // 16 for i in L]
    T13 = [(T10[i] - T9[i]) // 20 for i in L]
    T14 = [(T12[i] - T11[i]) // 21 for i in L]
    T15 = [(T13[i] - T12[i]) // 27 for i in L]
    r12 = [(T15[i] - T14[i]) // 22 for i in L]
    r10 = [T14[i] - 55*r12[i] for i in L]
    r8 = [T11[i] - 30*r10[i] - 627*r12[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] - 1408*r12[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] - 341*r12[i] 
          for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] - r12[i] for i in L]
    T16 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T17 = [((r[2][i] - r[-2][i]) // 4 - T16[i]) // 3 for i in L]
    T18 = [((r[3][i] - r[-3][i])//6 - T16[i])//8 for i in L]
    T19 = [((r[4][i] - r[-4][i])//8 - T16[i])//15 for i in L]
    T20 = [((r[5][i] - r[-5][i])//10 - T16[i])//24 for i in L]
    T21 = [((r[6][i] - r[-6][i])//12 - T16[i])//35 for i in L]
    T22 = [((r[7][i] - r0[i] - 49*r2[i] - 2401*r4[i] - 117649*r6[i] 
            - 5764801*r8[i] - 282475249*r10[i]
            - 13841287201*r12[i] - 678223072849*r14[i])//7 - T16[i])//48
            for i in L]
    T23 = [(T18[i] - T17[i]) // 5 for i in L]
    T24 = [(T19[i] - T18[i]) // 7 for i in L]
    T25 = [(T20[i] - T19[i]) // 9 for i in L]
    T26 = [(T21[i] - T20[i]) // 11 for i in L]
    T27 = [(T22[i] - T21[i]) // 13 for i in L]
    T28 = [(T24[i] - T23[i]) // 12 for i in L]
    T29 = [(T25[i] - T24[i]) // 16 for i in L]
    T30 = [(T26[i] - T25[i]) // 20 for i in L]
    T31 = [(T27[i] - T26[i]) // 24 for i in L]
    T32 = [(T29[i] - T28[i]) // 21 for i in L]
    T33 = [(T30[i] - T29[i]) // 27 for i in L]
    T34 = [(T31[i] - T30[i]) // 33 for i in L]
    T35 = [(T33[i] - T32[i]) // 32 for i in L]
    T36 = [(T34[i] - T33[i]) // 40 for i in L]
    r13 = [(T36[i] - T35[i]) // 45 for i in L]
    r11 = [T35[i] - 91*r13[i] for i in L]
    r9 = [T32[i] - 55*r11[i] - 2002*r13[i] for i in L]
    r7 = [T28[i] - 30*r9[i] - 627*r11[i] - 11440*r13[i] for i in L]
    r5 = [T23[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] - 13013*r13[i]
          for i in L]
    r3 = [T17[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i] - 1365*r13[i]
          for i in L]
    r1 = [T16[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] - r13[i]
          for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom9(f, g):
    n = 9
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, -6, 7, -7, 8, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]

    r16 = r['infinity']

    r14 = [(3003*(r[1][i] + r[-1][i])
     - 2002*(r[2][i] + r[-2][i])
     + 1001*(r[3][i] + r[-3][i])
     - 364*(r[4][i] + r[-4][i])
     + 91*(r[5][i] + r[-5][i])
     - 14*(r[6][i] + r[-6][i])
     + (r[7][i] + r[-7][i])
     - 3432*r0[i]
     - 12204960768000*r16[i]
    ) // 87178291200 for i in range(len(r0))]

    r12 = [(-792*(r[1][i] + r[-1][i])
     + 495*(r[2][i] + r[-2][i])
     - 220*(r[3][i] + r[-3][i])
     + 66*(r[4][i] + r[-4][i])
     - 12*(r[5][i] + r[-5][i])
     + (r[6][i] + r[-6][i])
     + 924*r0[i]
     - 43589145600*r14[i]
     - 2528170444800*r16[i]
    ) // 479001600 for i in range(len(r0))]

    r10 = [(210*(r[1][i] + r[-1][i])
     - 120*(r[2][i] + r[-2][i])
     + 45*(r[3][i] + r[-3][i])
     - 10*(r[4][i] + r[-4][i])
     + (r[5][i] + r[-5][i])
     - 252*r0[i]
     - 199584000*r12[i]
     - 7264857600*r14[i]
     - 223134912000*r16[i]
    ) // 3628800 for i in range(len(r0))]

    r8 = [(-56*(r[1][i] + r[-1][i])
     + 28*(r[2][i] + r[-2][i])
     - 8*(r[3][i] + r[-3][i])
     + (r[4][i] + r[-4][i])
     + 70*r0[i]
     - 1209600*r10[i]
     - 25280640*r12[i]
     - 461260800*r14[i]
     - 7904856960*r16[i]
    ) // 40320 for i in range(len(r0))]

    r6 = [(15*(r[1][i] + r[-1][i])
     - 6*(r[2][i] + r[-2][i])
     + (r[3][i] + r[-3][i])
     - 20*r0[i]
     - 10080*r8[i]
     - 105840*r10[i]
     - 1013760*r12[i]
     - 9369360*r14[i]
     - 85307040*r16[i]
    ) // 720 for i in range(len(r0))]

    r4 = [(-4*(r[1][i] + r[-1][i])
     + (r[2][i] + r[-2][i])
     + 6*r0[i]
     - 120*r6[i]
     - 504*r8[i]
     - 2040*r10[i]
     - 8184*r12[i]
     - 32760*r14[i]
     - 131064*r16[i]
    ) // 24 for i in range(len(r0))]

    r2 = [((r[1][i] + r[-1][i])
     - 2*r0[i]
     - 2*r4[i]
     - 2*r6[i]
     - 2*r8[i]
     - 2*r10[i]
     - 2*r12[i]
     - 2*r14[i]
     - 2*r16[i]
    ) // 2 for i in range(len(r0))]

    r15 = [(-1430*r[1][i] 
     + 2002*r[2][i]
     - 1638*r[3][i]
     + 910*r[4][i]
     - 350*r[5][i]
     + 90*r[6][i]
     - 14*r[7][i]
     + r[8][i]
     + 429*r0[i]
     - 264*r2[i]
     + 744*r4[i]
     - 5304*r6[i]
     + 81384*r8[i]
     - 2605944*r10[i]
     + 192387624*r12[i]
     - 55942352184*r14[i]
     - 20546119600536*r16[i]
    ) // 1307674368000 for i in range(len(r0))]

    r13 = [(429*r[1][i] 
     - 572*r[2][i]
     + 429*r[3][i]
     - 208*r[4][i]
     + 65*r[5][i]
     - 12*r[6][i]
     + r[7][i]
     - 132*r0[i]
     + 84*r2[i]
     - 252*r4[i]
     + 2004*r6[i]
     - 37212*r8[i]
     + 1710324*r10[i]
     - 325024572*r12[i]
     - 80789566956*r14[i]
     - 871782912000*r15[i]
     - 8422900930332*r16[i]
    ) // 6227020800 for i in range(len(r0))]

    r11 = [(-132*r[1][i] 
     + 165*r[2][i]
     - 110*r[3][i]
     + 44*r[4][i]
     - 10*r[5][i]
     + r[6][i]
     + 42*r0[i]
     - 28*r2[i]
     + 92*r4[i]
     - 868*r6[i]
     + 22652*r8[i]
     - 2620708*r10[i]
     - 415790788*r12[i]
     - 3632428800*r13[i]
     - 28616744548*r14[i]
     - 210680870400*r15[i]
     - 1479485236228*r16[i]
    ) // 39916800 for i in range(len(r0))]

    r9 = [(42*r[1][i] 
     - 48*r[2][i]
     + 27*r[3][i]
     - 8*r[4][i]
     + r[5][i]
     - 14*r0[i]
     + 10*r2[i]
     - 38*r4[i]
     + 490*r6[i]
     - 31238*r8[i]
     - 2922230*r10[i]
     - 19958400*r11[i]
     - 124075238*r12[i]
     - 726485760*r13[i]
     - 4084385750*r14[i]
     - 22313491200*r15[i]
     - 119387268038*r16[i]
    ) // 362880 for i in range(len(r0))]

    r7 = [(-14*r[1][i] 
     + 14*r[2][i]
     - 6*r[3][i]
     + r[4][i]
     + 5*r0[i]
     - 4*r2[i]
     + 20*r4[i]
     - 604*r6[i]
     - 29740*r8[i]
     - 151200*r9[i]
     - 708604*r10[i]
     - 3160080*r11[i]
     - 13645900*r12[i]
     - 57657600*r13[i]
     - 239967004*r14[i]
     - 988107120*r15[i]
     - 4037604460*r16[i]
    ) // 5040 for i in range(len(r0))]

    r5 = [(5*r[1][i] 
     - 4*r[2][i]
     + r[3][i]
     - 2*r0[i]
     + 2*r2[i]
     - 22*r4[i]
     - 478*r6[i]
     - 1680*r7[i]
     - 5542*r8[i]
     - 17640*r9[i]
     - 54958*r10[i]
     - 168960*r11[i]
     - 515062*r12[i]
     - 1561560*r13[i]
     - 4717438*r14[i]
     - 14217840*r15[i]
     - 42784582*r16[i]
    ) // 120 for i in range(len(r0))]

    r3 = [(-2*r[1][i] 
     + r[2][i]
     + 1*r0[i]
     - 2*r2[i]
     - 14*r4[i]
     - 30*r5[i]
     - 62*r6[i]
     - 126*r7[i]
     - 254*r8[i]
     - 510*r9[i]
     - 1022*r10[i]
     - 2046*r11[i]
     - 4094*r12[i]
     - 8190*r13[i]
     - 16382*r14[i]
     - 32766*r15[i]
     - 65534*r16[i]
    ) // 6 for i in range(len(r0))]

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
    - r11[i]
    - r12[i]
    - r13[i]
    - r14[i]
    - r15[i]
    - r16[i]
     for i in range(len(r0))]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14,
               r15, r16)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]
# =======================================
#
#          Compound Algorithms
#
# =======================================

def toom2_toom2(f, g):
    (f0, f1) = split(f, 2)
    (g0, g1) = split(g, 2)
    
    r0 = toom2(f0, g0)
    r2 = toom2(f1, g1)
    temp_prod = toom2(poly_add(f0, f1), poly_add(g0, g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1] 

def toom2_toom2_toom2(f, g):
    (f0, f1) = split(f, 2)
    (g0, g1) = split(g, 2)
    
    r0 = toom2_toom2(f0, g0)
    r2 = toom2_toom2(f1, g1)
    temp_prod = toom2_toom2(poly_add(f0, f1), poly_add(g0, g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1]

def toom2_toom2_toom2_toom2(f, g):
    (f0, f1) = split(f, 2)
    (g0, g1) = split(g, 2)
    
    r0 = toom2_toom2_toom2(f0, g0)
    r2 = toom2_toom2_toom2(f1, g1)
    temp_prod = toom2_toom2_toom2(poly_add(f0, f1), poly_add(g0, g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1] 

def toom2_toom4_toom3(f, g):
    (f0, f1) = split(f, 2)
    (g0, g1) = split(g, 2)
    
    r0 = toom4_toom3(f0, g0)
    r2 = toom4_toom3(f1, g1)
    temp_prod = toom4_toom3(poly_add(f0, f1), poly_add(g0, g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1] 

def toom3_toom2(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
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

def toom3_toom2_toom2(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2_toom2(f_eval[i], g_eval[i])
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

def toom3_toom3(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
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

def toom3_toom3_toom2(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3_toom2(f_eval[i], g_eval[i])
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

def toom3_toom2_toom2_toom2(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2_toom2_toom2(f_eval[i], g_eval[i])
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

def toom3_toom3_toom3(f, g):
    n = 3
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3_toom3(f_eval[i], g_eval[i])
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


def toom4_toom2(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
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

def toom4_toom3(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
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


def toom4_toom4(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom4(f_eval[i], g_eval[i])
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


def toom4_toom2_toom2(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2_toom2(f_eval[i], g_eval[i])
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


def toom4_toom3_toom2(f, g):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3_toom2(f_eval[i], g_eval[i])
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

def toom5_toom2(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom5_toom3(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom5_toom4(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom4(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom5_toom5(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom5(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom5_toom2_toom2(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2_toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom5_toom3_toom2(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3_toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r8 = r['infinity']
    L = range(len(r0)) # to save time on for-loops
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r8[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i]) // 2 - r0[i] - 256*r8[i]) // 4 - T1[i]) // 3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i]) // 2 - r0[i] - 6561*r8[i]) // 9 - T1[i]) // 8
          for i in L]
    r6 = [(T3[i] - T2[i]) // 5 for i in L]
    r4 = [T2[i] - 5*r6[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] for i in L]
    T4 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T5 = [((r[2][i] - r[-2][i]) // 4 - T4[i]) // 3 for i in L]
    T6 = [((r[3][i] - r[-3][i]) // 6 - T4[i]) // 8 for i in L]
    T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])
          // 4 - T4[i]) // 15 for i in L]
    T8 = [(T6[i] - T5[i]) // 5 for i in L]
    T9 = [(T7[i] - T6[i]) // 7 for i in L]
    r7 = [(T9[i] - T8[i]) // 12 for i in L]
    r5 = [T8[i] - 14*r7[i] for i in L]
    r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in L]
    r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom6_toom2(f, g):
    n = 6
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r10 = r['infinity']
    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i]
          for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 1024*r10[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 59049*r10[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 1048576*r10[i])//16 - T1[i])//15
          for i in L]
    T5 = [(T3[i] - T2[i])//5 for i in L]
    T6 = [(T4[i] - T3[i])//7 for i in L]
    r8 = [(T6[i] - T5[i])//12 for i in L]
    r6 = [T5[i] - 14*r8[i] for i in L]
    r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in L]
    
    T7 = [(r[1][i] - r[-1][i])//2 for i in L]
    T8 = [((r[2][i] - r[-2][i])//4 - T7[i])//3 for i in L]
    T9 = [((r[3][i] - r[-3][i])//6 - T7[i])//8 for i in L]
    T10 = [((r[4][i] - r[-4][i])//8 - T7[i])//15 for i in L]
    T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] - 15625*r6[i] 
    - 390625*r8[i] - 9765625*r10[i])//5 - T7[i])//24 for i in L]
    T12 = [(T9[i] - T8[i])//5 for i in L]
    T13 = [(T10[i] - T9[i])//7 for i in L]
    T14 = [(T11[i] - T10[i])//9 for i in L]
    T15 = [(T13[i] - T12[i])//12 for i in L]
    T16 = [(T14[i] - T13[i])//16 for i in L]
    r9 = [(T16[i] - T15[i])//21 for i in L]
    r7 = [T15[i] - 30*r9[i] for i in L]
    r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in L]
    r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in L]
    r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom6_toom4(f, g):
    n = 6
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom4(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r10 = r['infinity']
    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i]
          for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 1024*r10[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 59049*r10[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 1048576*r10[i])//16 - T1[i])//15
          for i in L]
    T5 = [(T3[i] - T2[i])//5 for i in L]
    T6 = [(T4[i] - T3[i])//7 for i in L]
    r8 = [(T6[i] - T5[i])//12 for i in L]
    r6 = [T5[i] - 14*r8[i] for i in L]
    r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in L]
    
    T7 = [(r[1][i] - r[-1][i])//2 for i in L]
    T8 = [((r[2][i] - r[-2][i])//4 - T7[i])//3 for i in L]
    T9 = [((r[3][i] - r[-3][i])//6 - T7[i])//8 for i in L]
    T10 = [((r[4][i] - r[-4][i])//8 - T7[i])//15 for i in L]
    T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] - 15625*r6[i] 
    - 390625*r8[i] - 9765625*r10[i])//5 - T7[i])//24 for i in L]
    T12 = [(T9[i] - T8[i])//5 for i in L]
    T13 = [(T10[i] - T9[i])//7 for i in L]
    T14 = [(T11[i] - T10[i])//9 for i in L]
    T15 = [(T13[i] - T12[i])//12 for i in L]
    T16 = [(T14[i] - T13[i])//16 for i in L]
    r9 = [(T16[i] - T15[i])//21 for i in L]
    r7 = [T15[i] - 30*r9[i] for i in L]
    r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in L]
    r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in L]
    r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom6_toom5(f, g):
    n = 6
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom5(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r10 = r['infinity']
    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i]
          for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 1024*r10[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 59049*r10[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 1048576*r10[i])//16 - T1[i])//15
          for i in L]
    T5 = [(T3[i] - T2[i])//5 for i in L]
    T6 = [(T4[i] - T3[i])//7 for i in L]
    r8 = [(T6[i] - T5[i])//12 for i in L]
    r6 = [T5[i] - 14*r8[i] for i in L]
    r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in L]
    
    T7 = [(r[1][i] - r[-1][i])//2 for i in L]
    T8 = [((r[2][i] - r[-2][i])//4 - T7[i])//3 for i in L]
    T9 = [((r[3][i] - r[-3][i])//6 - T7[i])//8 for i in L]
    T10 = [((r[4][i] - r[-4][i])//8 - T7[i])//15 for i in L]
    T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] - 15625*r6[i] 
    - 390625*r8[i] - 9765625*r10[i])//5 - T7[i])//24 for i in L]
    T12 = [(T9[i] - T8[i])//5 for i in L]
    T13 = [(T10[i] - T9[i])//7 for i in L]
    T14 = [(T11[i] - T10[i])//9 for i in L]
    T15 = [(T13[i] - T12[i])//12 for i in L]
    T16 = [(T14[i] - T13[i])//16 for i in L]
    r9 = [(T16[i] - T15[i])//21 for i in L]
    r7 = [T15[i] - 30*r9[i] for i in L]
    r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in L]
    r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in L]
    r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom7_toom2(f, g):
    n = 7
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
    L = range(len(r0))
    r12 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r12[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 4096*r12[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 531441*r12[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 16777216*r12[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 244140625*r12[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(T3[i] - T2[i]) // 5 for i in L]
    T7 = [(T4[i] - T3[i]) // 7 for i in L]
    T8 = [(T5[i] - T4[i]) // 9 for i in L]
    T9 = [(T7[i] - T6[i]) // 12 for i in L]
    T10 = [(T8[i] - T7[i]) // 16 for i in L]
    r10 = [(T10[i] - T9[i]) // 21 for i in L]
    r8 = [T9[i] - 30*r10[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] for i in L]
    T11 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T12 = [((r[2][i] - r[-2][i]) // 4 - T11[i]) // 3 for i in L]
    T13 = [((r[3][i] - r[-3][i])//6 - T11[i])//8 for i in L]
    T14 = [((r[4][i] - r[-4][i])//8 - T11[i])//15 for i in L]
    T15 = [((r[5][i] - r[-5][i])//10 - T11[i])//24 for i in L]
    T16 = [((r[6][i] - r0[i] - 36*r2[i] - 1296*r4[i] - 46656*r6[i] 
           - 1679616*r8[i] - 60466176*r10[i] - 2176782336*r12[i])//6 
           - T11[i])//35 for i in L]
    T17 = [(T13[i] - T12[i]) // 5 for i in L]
    T18 = [(T14[i] - T13[i]) // 7 for i in L]
    T19 = [(T15[i] - T14[i]) // 9 for i in L]
    T20 = [(T16[i] - T15[i]) // 11 for i in L]
    T21 = [(T18[i] - T17[i]) // 12 for i in L]
    T22 = [(T19[i] - T18[i]) // 16 for i in L]
    T23 = [(T20[i] - T19[i]) // 20 for i in L]
    T24 = [(T22[i] - T21[i]) // 21 for i in L]
    T25 = [(T23[i] - T22[i]) // 27 for i in L]
    r11 = [(T25[i] - T24[i]) // 32 for i in L]
    r9 = [T24[i] - 55*r11[i] for i in L]
    r7 = [T21[i] - 30*r9[i] - 627*r11[i] for i in L]
    r5 = [T17[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] for i in L]
    r3 = [T12[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i]
          for i in L]
    r1 = [T11[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom7_toom3(f, g):
    n = 7
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
    L = range(len(r0))
    r12 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r12[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 4096*r12[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 531441*r12[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 16777216*r12[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 244140625*r12[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(T3[i] - T2[i]) // 5 for i in L]
    T7 = [(T4[i] - T3[i]) // 7 for i in L]
    T8 = [(T5[i] - T4[i]) // 9 for i in L]
    T9 = [(T7[i] - T6[i]) // 12 for i in L]
    T10 = [(T8[i] - T7[i]) // 16 for i in L]
    r10 = [(T10[i] - T9[i]) // 21 for i in L]
    r8 = [T9[i] - 30*r10[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] for i in L]
    T11 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T12 = [((r[2][i] - r[-2][i]) // 4 - T11[i]) // 3 for i in L]
    T13 = [((r[3][i] - r[-3][i])//6 - T11[i])//8 for i in L]
    T14 = [((r[4][i] - r[-4][i])//8 - T11[i])//15 for i in L]
    T15 = [((r[5][i] - r[-5][i])//10 - T11[i])//24 for i in L]
    T16 = [((r[6][i] - r0[i] - 36*r2[i] - 1296*r4[i] - 46656*r6[i] 
           - 1679616*r8[i] - 60466176*r10[i] - 2176782336*r12[i])//6 
           - T11[i])//35 for i in L]
    T17 = [(T13[i] - T12[i]) // 5 for i in L]
    T18 = [(T14[i] - T13[i]) // 7 for i in L]
    T19 = [(T15[i] - T14[i]) // 9 for i in L]
    T20 = [(T16[i] - T15[i]) // 11 for i in L]
    T21 = [(T18[i] - T17[i]) // 12 for i in L]
    T22 = [(T19[i] - T18[i]) // 16 for i in L]
    T23 = [(T20[i] - T19[i]) // 20 for i in L]
    T24 = [(T22[i] - T21[i]) // 21 for i in L]
    T25 = [(T23[i] - T22[i]) // 27 for i in L]
    r11 = [(T25[i] - T24[i]) // 32 for i in L]
    r9 = [T24[i] - 55*r11[i] for i in L]
    r7 = [T21[i] - 30*r9[i] - 627*r11[i] for i in L]
    r5 = [T17[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] for i in L]
    r3 = [T12[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i]
          for i in L]
    r1 = [T11[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]



def toom7_toom2_toom2(f, g):
    n = 7
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i] : toom2_toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
    L = range(len(r0))
    r12 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r12[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 4096*r12[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 531441*r12[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 16777216*r12[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 244140625*r12[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(T3[i] - T2[i]) // 5 for i in L]
    T7 = [(T4[i] - T3[i]) // 7 for i in L]
    T8 = [(T5[i] - T4[i]) // 9 for i in L]
    T9 = [(T7[i] - T6[i]) // 12 for i in L]
    T10 = [(T8[i] - T7[i]) // 16 for i in L]
    r10 = [(T10[i] - T9[i]) // 21 for i in L]
    r8 = [T9[i] - 30*r10[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] for i in L]
    T11 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T12 = [((r[2][i] - r[-2][i]) // 4 - T11[i]) // 3 for i in L]
    T13 = [((r[3][i] - r[-3][i])//6 - T11[i])//8 for i in L]
    T14 = [((r[4][i] - r[-4][i])//8 - T11[i])//15 for i in L]
    T15 = [((r[5][i] - r[-5][i])//10 - T11[i])//24 for i in L]
    T16 = [((r[6][i] - r0[i] - 36*r2[i] - 1296*r4[i] - 46656*r6[i] 
           - 1679616*r8[i] - 60466176*r10[i] - 2176782336*r12[i])//6 
           - T11[i])//35 for i in L]
    T17 = [(T13[i] - T12[i]) // 5 for i in L]
    T18 = [(T14[i] - T13[i]) // 7 for i in L]
    T19 = [(T15[i] - T14[i]) // 9 for i in L]
    T20 = [(T16[i] - T15[i]) // 11 for i in L]
    T21 = [(T18[i] - T17[i]) // 12 for i in L]
    T22 = [(T19[i] - T18[i]) // 16 for i in L]
    T23 = [(T20[i] - T19[i]) // 20 for i in L]
    T24 = [(T22[i] - T21[i]) // 21 for i in L]
    T25 = [(T23[i] - T22[i]) // 27 for i in L]
    r11 = [(T25[i] - T24[i]) // 32 for i in L]
    r9 = [T24[i] - 55*r11[i] for i in L]
    r7 = [T21[i] - 30*r9[i] - 627*r11[i] for i in L]
    r5 = [T17[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] for i in L]
    r3 = [T12[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i]
          for i in L]
    r1 = [T11[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom7_toom4(f, g):
    n = 7
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i] : toom4(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]
    L = range(len(r0))
    r12 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r12[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 4096*r12[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 531441*r12[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 16777216*r12[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 244140625*r12[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(T3[i] - T2[i]) // 5 for i in L]
    T7 = [(T4[i] - T3[i]) // 7 for i in L]
    T8 = [(T5[i] - T4[i]) // 9 for i in L]
    T9 = [(T7[i] - T6[i]) // 12 for i in L]
    T10 = [(T8[i] - T7[i]) // 16 for i in L]
    r10 = [(T10[i] - T9[i]) // 21 for i in L]
    r8 = [T9[i] - 30*r10[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] for i in L]
    T11 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T12 = [((r[2][i] - r[-2][i]) // 4 - T11[i]) // 3 for i in L]
    T13 = [((r[3][i] - r[-3][i])//6 - T11[i])//8 for i in L]
    T14 = [((r[4][i] - r[-4][i])//8 - T11[i])//15 for i in L]
    T15 = [((r[5][i] - r[-5][i])//10 - T11[i])//24 for i in L]
    T16 = [((r[6][i] - r0[i] - 36*r2[i] - 1296*r4[i] - 46656*r6[i] 
           - 1679616*r8[i] - 60466176*r10[i] - 2176782336*r12[i])//6 
           - T11[i])//35 for i in L]
    T17 = [(T13[i] - T12[i]) // 5 for i in L]
    T18 = [(T14[i] - T13[i]) // 7 for i in L]
    T19 = [(T15[i] - T14[i]) // 9 for i in L]
    T20 = [(T16[i] - T15[i]) // 11 for i in L]
    T21 = [(T18[i] - T17[i]) // 12 for i in L]
    T22 = [(T19[i] - T18[i]) // 16 for i in L]
    T23 = [(T20[i] - T19[i]) // 20 for i in L]
    T24 = [(T22[i] - T21[i]) // 21 for i in L]
    T25 = [(T23[i] - T22[i]) // 27 for i in L]
    r11 = [(T25[i] - T24[i]) // 32 for i in L]
    r9 = [T24[i] - 55*r11[i] for i in L]
    r7 = [T21[i] - 30*r9[i] - 627*r11[i] for i in L]
    r5 = [T17[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] for i in L]
    r3 = [T12[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i]
          for i in L]
    r1 = [T11[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]



def toom8_toom2(f, g):
    n = 8
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, -6, 7, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom2(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r14 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r14[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 16384*r14[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 4782969*r14[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 268435456*r14[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 6103515625*r14[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(((r[6][i] + r[-6][i])//2 - r0[i] - 78364164096*r14[i])//36 
           - T1[i])//35 for i in L]
    T7 = [(T3[i] - T2[i]) // 5 for i in L]
    T8 = [(T4[i] - T3[i]) // 7 for i in L]
    T9 = [(T5[i] - T4[i]) // 9 for i in L]
    T10 = [(T6[i] - T5[i]) // 11 for i in L]
    T11 = [(T8[i] - T7[i]) // 12 for i in L]
    T12 = [(T9[i] - T8[i]) // 16 for i in L]
    T13 = [(T10[i] - T9[i]) // 20 for i in L]
    T14 = [(T12[i] - T11[i]) // 21 for i in L]
    T15 = [(T13[i] - T12[i]) // 27 for i in L]
    r12 = [(T15[i] - T14[i]) // 22 for i in L]
    r10 = [T14[i] - 55*r12[i] for i in L]
    r8 = [T11[i] - 30*r10[i] - 627*r12[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] - 1408*r12[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] - 341*r12[i] 
          for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] - r12[i] for i in L]
    T16 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T17 = [((r[2][i] - r[-2][i]) // 4 - T16[i]) // 3 for i in L]
    T18 = [((r[3][i] - r[-3][i])//6 - T16[i])//8 for i in L]
    T19 = [((r[4][i] - r[-4][i])//8 - T16[i])//15 for i in L]
    T20 = [((r[5][i] - r[-5][i])//10 - T16[i])//24 for i in L]
    T21 = [((r[6][i] - r[-6][i])//12 - T16[i])//35 for i in L]
    T22 = [((r[7][i] - r0[i] - 49*r2[i] - 2401*r4[i] - 117649*r6[i] 
            - 5764801*r8[i] - 282475249*r10[i]
            - 13841287201*r12[i] - 678223072849*r14[i])//7 - T16[i])//48
            for i in L]
    T23 = [(T18[i] - T17[i]) // 5 for i in L]
    T24 = [(T19[i] - T18[i]) // 7 for i in L]
    T25 = [(T20[i] - T19[i]) // 9 for i in L]
    T26 = [(T21[i] - T20[i]) // 11 for i in L]
    T27 = [(T22[i] - T21[i]) // 13 for i in L]
    T28 = [(T24[i] - T23[i]) // 12 for i in L]
    T29 = [(T25[i] - T24[i]) // 16 for i in L]
    T30 = [(T26[i] - T25[i]) // 20 for i in L]
    T31 = [(T27[i] - T26[i]) // 24 for i in L]
    T32 = [(T29[i] - T28[i]) // 21 for i in L]
    T33 = [(T30[i] - T29[i]) // 27 for i in L]
    T34 = [(T31[i] - T30[i]) // 33 for i in L]
    T35 = [(T33[i] - T32[i]) // 32 for i in L]
    T36 = [(T34[i] - T33[i]) // 40 for i in L]
    r13 = [(T36[i] - T35[i]) // 45 for i in L]
    r11 = [T35[i] - 91*r13[i] for i in L]
    r9 = [T32[i] - 55*r11[i] - 2002*r13[i] for i in L]
    r7 = [T28[i] - 30*r9[i] - 627*r11[i] - 11440*r13[i] for i in L]
    r5 = [T23[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] - 13013*r13[i]
          for i in L]
    r3 = [T17[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i] - 1365*r13[i]
          for i in L]
    r1 = [T16[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] - r13[i]
          for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def toom8_toom3(f, g):
    n = 8
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, -6, 7, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    L = range(len(r0))
    r14 = r['infinity']
    T1 = [(r[1][i] + r[-1][i]) // 2 - r0[i] - r14[i] for i in L]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 16384*r14[i])//4 - T1[i])//3
          for i in L]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 4782969*r14[i])//9 - T1[i])//8
          for i in L]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 268435456*r14[i])//16 
           - T1[i])//15 for i in L]
    T5 = [(((r[5][i] + r[-5][i])//2 - r0[i] - 6103515625*r14[i])//25 
           - T1[i])//24 for i in L]
    T6 = [(((r[6][i] + r[-6][i])//2 - r0[i] - 78364164096*r14[i])//36 
           - T1[i])//35 for i in L]
    T7 = [(T3[i] - T2[i]) // 5 for i in L]
    T8 = [(T4[i] - T3[i]) // 7 for i in L]
    T9 = [(T5[i] - T4[i]) // 9 for i in L]
    T10 = [(T6[i] - T5[i]) // 11 for i in L]
    T11 = [(T8[i] - T7[i]) // 12 for i in L]
    T12 = [(T9[i] - T8[i]) // 16 for i in L]
    T13 = [(T10[i] - T9[i]) // 20 for i in L]
    T14 = [(T12[i] - T11[i]) // 21 for i in L]
    T15 = [(T13[i] - T12[i]) // 27 for i in L]
    r12 = [(T15[i] - T14[i]) // 22 for i in L]
    r10 = [T14[i] - 55*r12[i] for i in L]
    r8 = [T11[i] - 30*r10[i] - 627*r12[i] for i in L]
    r6 = [T7[i] - 14*r8[i] - 147*r10[i] - 1408*r12[i] for i in L]
    r4 = [T2[i] - 5*r8[i] - 21*r8[i] - 85*r10[i] - 341*r12[i] 
          for i in L]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] - r10[i] - r12[i] for i in L]
    T16 = [(r[1][i] - r[-1][i]) // 2 for i in L]
    T17 = [((r[2][i] - r[-2][i]) // 4 - T16[i]) // 3 for i in L]
    T18 = [((r[3][i] - r[-3][i])//6 - T16[i])//8 for i in L]
    T19 = [((r[4][i] - r[-4][i])//8 - T16[i])//15 for i in L]
    T20 = [((r[5][i] - r[-5][i])//10 - T16[i])//24 for i in L]
    T21 = [((r[6][i] - r[-6][i])//12 - T16[i])//35 for i in L]
    T22 = [((r[7][i] - r0[i] - 49*r2[i] - 2401*r4[i] - 117649*r6[i] 
            - 5764801*r8[i] - 282475249*r10[i]
            - 13841287201*r12[i] - 678223072849*r14[i])//7 - T16[i])//48
            for i in L]
    T23 = [(T18[i] - T17[i]) // 5 for i in L]
    T24 = [(T19[i] - T18[i]) // 7 for i in L]
    T25 = [(T20[i] - T19[i]) // 9 for i in L]
    T26 = [(T21[i] - T20[i]) // 11 for i in L]
    T27 = [(T22[i] - T21[i]) // 13 for i in L]
    T28 = [(T24[i] - T23[i]) // 12 for i in L]
    T29 = [(T25[i] - T24[i]) // 16 for i in L]
    T30 = [(T26[i] - T25[i]) // 20 for i in L]
    T31 = [(T27[i] - T26[i]) // 24 for i in L]
    T32 = [(T29[i] - T28[i]) // 21 for i in L]
    T33 = [(T30[i] - T29[i]) // 27 for i in L]
    T34 = [(T31[i] - T30[i]) // 33 for i in L]
    T35 = [(T33[i] - T32[i]) // 32 for i in L]
    T36 = [(T34[i] - T33[i]) // 40 for i in L]
    r13 = [(T36[i] - T35[i]) // 45 for i in L]
    r11 = [T35[i] - 91*r13[i] for i in L]
    r9 = [T32[i] - 55*r11[i] - 2002*r13[i] for i in L]
    r7 = [T28[i] - 30*r9[i] - 627*r11[i] - 11440*r13[i] for i in L]
    r5 = [T23[i] - 14*r7[i] - 147*r9[i] - 1408*r11[i] - 13013*r13[i]
          for i in L]
    r3 = [T17[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] - 341*r11[i] - 1365*r13[i]
          for i in L]
    r1 = [T16[i] - r3[i] - r5[i] - r7[i] - r9[i] - r11[i] - r13[i]
          for i in L]
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def toom9_toom3(f, g):
    n = 9
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5,
                  -5, 6, -6, 7, -7, 8, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)

    r = {eval_list[i] : toom3(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}

    r0 = r[0]

    r16 = r['infinity']

    r14 = [(3003*(r[1][i] + r[-1][i])
     - 2002*(r[2][i] + r[-2][i])
     + 1001*(r[3][i] + r[-3][i])
     - 364*(r[4][i] + r[-4][i])
     + 91*(r[5][i] + r[-5][i])
     - 14*(r[6][i] + r[-6][i])
     + (r[7][i] + r[-7][i])
     - 3432*r0[i]
     - 12204960768000*r16[i]
    ) // 87178291200 for i in range(len(r0))]

    r12 = [(-792*(r[1][i] + r[-1][i])
     + 495*(r[2][i] + r[-2][i])
     - 220*(r[3][i] + r[-3][i])
     + 66*(r[4][i] + r[-4][i])
     - 12*(r[5][i] + r[-5][i])
     + (r[6][i] + r[-6][i])
     + 924*r0[i]
     - 43589145600*r14[i]
     - 2528170444800*r16[i]
    ) // 479001600 for i in range(len(r0))]

    r10 = [(210*(r[1][i] + r[-1][i])
     - 120*(r[2][i] + r[-2][i])
     + 45*(r[3][i] + r[-3][i])
     - 10*(r[4][i] + r[-4][i])
     + (r[5][i] + r[-5][i])
     - 252*r0[i]
     - 199584000*r12[i]
     - 7264857600*r14[i]
     - 223134912000*r16[i]
    ) // 3628800 for i in range(len(r0))]

    r8 = [(-56*(r[1][i] + r[-1][i])
     + 28*(r[2][i] + r[-2][i])
     - 8*(r[3][i] + r[-3][i])
     + (r[4][i] + r[-4][i])
     + 70*r0[i]
     - 1209600*r10[i]
     - 25280640*r12[i]
     - 461260800*r14[i]
     - 7904856960*r16[i]
    ) // 40320 for i in range(len(r0))]

    r6 = [(15*(r[1][i] + r[-1][i])
     - 6*(r[2][i] + r[-2][i])
     + (r[3][i] + r[-3][i])
     - 20*r0[i]
     - 10080*r8[i]
     - 105840*r10[i]
     - 1013760*r12[i]
     - 9369360*r14[i]
     - 85307040*r16[i]
    ) // 720 for i in range(len(r0))]

    r4 = [(-4*(r[1][i] + r[-1][i])
     + (r[2][i] + r[-2][i])
     + 6*r0[i]
     - 120*r6[i]
     - 504*r8[i]
     - 2040*r10[i]
     - 8184*r12[i]
     - 32760*r14[i]
     - 131064*r16[i]
    ) // 24 for i in range(len(r0))]

    r2 = [((r[1][i] + r[-1][i])
     - 2*r0[i]
     - 2*r4[i]
     - 2*r6[i]
     - 2*r8[i]
     - 2*r10[i]
     - 2*r12[i]
     - 2*r14[i]
     - 2*r16[i]
    ) // 2 for i in range(len(r0))]

    r15 = [(-1430*r[1][i] 
     + 2002*r[2][i]
     - 1638*r[3][i]
     + 910*r[4][i]
     - 350*r[5][i]
     + 90*r[6][i]
     - 14*r[7][i]
     + r[8][i]
     + 429*r0[i]
     - 264*r2[i]
     + 744*r4[i]
     - 5304*r6[i]
     + 81384*r8[i]
     - 2605944*r10[i]
     + 192387624*r12[i]
     - 55942352184*r14[i]
     - 20546119600536*r16[i]
    ) // 1307674368000 for i in range(len(r0))]

    r13 = [(429*r[1][i] 
     - 572*r[2][i]
     + 429*r[3][i]
     - 208*r[4][i]
     + 65*r[5][i]
     - 12*r[6][i]
     + r[7][i]
     - 132*r0[i]
     + 84*r2[i]
     - 252*r4[i]
     + 2004*r6[i]
     - 37212*r8[i]
     + 1710324*r10[i]
     - 325024572*r12[i]
     - 80789566956*r14[i]
     - 871782912000*r15[i]
     - 8422900930332*r16[i]
    ) // 6227020800 for i in range(len(r0))]

    r11 = [(-132*r[1][i] 
     + 165*r[2][i]
     - 110*r[3][i]
     + 44*r[4][i]
     - 10*r[5][i]
     + r[6][i]
     + 42*r0[i]
     - 28*r2[i]
     + 92*r4[i]
     - 868*r6[i]
     + 22652*r8[i]
     - 2620708*r10[i]
     - 415790788*r12[i]
     - 3632428800*r13[i]
     - 28616744548*r14[i]
     - 210680870400*r15[i]
     - 1479485236228*r16[i]
    ) // 39916800 for i in range(len(r0))]

    r9 = [(42*r[1][i] 
     - 48*r[2][i]
     + 27*r[3][i]
     - 8*r[4][i]
     + r[5][i]
     - 14*r0[i]
     + 10*r2[i]
     - 38*r4[i]
     + 490*r6[i]
     - 31238*r8[i]
     - 2922230*r10[i]
     - 19958400*r11[i]
     - 124075238*r12[i]
     - 726485760*r13[i]
     - 4084385750*r14[i]
     - 22313491200*r15[i]
     - 119387268038*r16[i]
    ) // 362880 for i in range(len(r0))]

    r7 = [(-14*r[1][i] 
     + 14*r[2][i]
     - 6*r[3][i]
     + r[4][i]
     + 5*r0[i]
     - 4*r2[i]
     + 20*r4[i]
     - 604*r6[i]
     - 29740*r8[i]
     - 151200*r9[i]
     - 708604*r10[i]
     - 3160080*r11[i]
     - 13645900*r12[i]
     - 57657600*r13[i]
     - 239967004*r14[i]
     - 988107120*r15[i]
     - 4037604460*r16[i]
    ) // 5040 for i in range(len(r0))]

    r5 = [(5*r[1][i] 
     - 4*r[2][i]
     + r[3][i]
     - 2*r0[i]
     + 2*r2[i]
     - 22*r4[i]
     - 478*r6[i]
     - 1680*r7[i]
     - 5542*r8[i]
     - 17640*r9[i]
     - 54958*r10[i]
     - 168960*r11[i]
     - 515062*r12[i]
     - 1561560*r13[i]
     - 4717438*r14[i]
     - 14217840*r15[i]
     - 42784582*r16[i]
    ) // 120 for i in range(len(r0))]

    r3 = [(-2*r[1][i] 
     + r[2][i]
     + 1*r0[i]
     - 2*r2[i]
     - 14*r4[i]
     - 30*r5[i]
     - 62*r6[i]
     - 126*r7[i]
     - 254*r8[i]
     - 510*r9[i]
     - 1022*r10[i]
     - 2046*r11[i]
     - 4094*r12[i]
     - 8190*r13[i]
     - 16382*r14[i]
     - 32766*r15[i]
     - 65534*r16[i]
    ) // 6 for i in range(len(r0))]

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
    - r11[i]
    - r12[i]
    - r13[i]
    - r14[i]
    - r15[i]
    - r16[i]
     for i in range(len(r0))]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14,
               r15, r16)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def single_threshold(algorithm, 
                     num_trials=500, 
                     min_deg=5, 
                     max_deg=50, 
                     filename=None):
    """ This is for finding the threshold of a single Toom algorithm. Input
        an integer from 2 to 8?"""
    
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    # dictionaries keeping track of sums of times. We will divide them later
    # to get the averages
    tn_times = {degree:0 for degree in degrees}
    sb_times = {degree:0 for degree in degrees}
    
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
                
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            # time the schoolbook
            start_time = time.time()
            schoolbook(f, g)
            sb_times[degree] += time.time() - start_time
            
            # time the other algorithm
            if algorithm == 2:
                start_time = time.time()
                toom2(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 3:
                start_time = time.time()
                toom3(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 4:
                start_time = time.time()
                toom4(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 5:
                start_time = time.time()
                toom5(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 6:
                start_time = time.time()
                toom6(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 7:
                start_time = time.time()
                toom7(f, g)
                tn_times[degree] += time.time() - start_time
            elif algorithm == 8:
                start_time = time.time()
                toom8(f, g)
                tn_times[degree] += time.time() - start_time
            
        # progress bar
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")

    sb_times = [sb_times[degree]/num_trials for degree in degrees]
    tn_times = [tn_times[degree]/num_trials for degree in degrees]
    
    # names
    name = "Toom-{}".format(algorithm)
    
    plt.plot(degrees, sb_times, lw=3, label='Schoolbook')
    plt.plot(degrees, tn_times, lw=3, label=name)
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Time (seconds)", size=15)
    plt.title("Average Time for Polynomial Multiplication", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return


def plot_algorithms(algorithms, num_trials=100, min_deg=5, max_deg=50,
                    filename=None):
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    # dictionaries keeping track of sums of times. We will divide them later
    # to get the averages
    alg2time = {}
    for algorithm in algorithms:
        alg2time[algorithm] = {degree:0 for degree in degrees}
        
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
                
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            if 'sb' in algorithms:
                start_time = time.time()
                schoolbook(f, g)
                alg2time['sb'][degree] += time.time() - start_time
            if 2 in algorithms:
                start_time = time.time()
                toom2(f, g)
                alg2time[2][degree] += time.time() - start_time
            if 3 in algorithms:
                start_time = time.time()
                toom3(f, g)
                alg2time[3][degree] += time.time() - start_time
            if 4 in algorithms:
                start_time = time.time()
                toom4(f, g)
                alg2time[4][degree] += time.time() - start_time
            if 5 in algorithms:
                start_time = time.time()
                toom5(f, g)
                alg2time[5][degree] += time.time() - start_time
            if 6 in algorithms:
                start_time = time.time()
                toom6(f, g)
                alg2time[6][degree] += time.time() - start_time
            if 22 in algorithms:
                start_time = time.time()
                toom2_toom2(f, g)
                alg2time[22][degree] += time.time() - start_time
            if 32 in algorithms:
                start_time = time.time()
                toom3_toom2(f, g)
                alg2time[32][degree] += time.time() - start_time
            if 33 in algorithms:
                start_time = time.time()
                toom3_toom3(f, g)
                alg2time[33][degree] += time.time() - start_time
            if 222 in algorithms:
                start_time = time.time()
                toom2_toom2_toom2(f, g)
                alg2time[222][degree] += time.time() - start_time
            if 2222 in algorithms:
                start_time = time.time()
                toom2_toom2_toom2_toom2(f, g)
                alg2time[2222][degree] += time.time() - start_time
            if 42 in algorithms:
                start_time = time.time()
                toom4_toom2(f, g)
                alg2time[42][degree] += time.time() - start_time
            if 43 in algorithms:
                start_time = time.time()
                toom4_toom3(f, g)
            if 44 in algorithms:
                start_time = time.time()
                toom4_toom4(f, g)
                alg2time[44][degree] += time.time() - start_time
            if 422 in algorithms:
                start_time = time.time()
                toom4_toom2_toom2(f, g)
                alg2time[422][degree] += time.time() - start_time
            if 432 in algorithms:
                start_time = time.time()
                toom4_toom3_toom2(f, g)
                alg2time[432][degree] += time.time() - start_time
            if 243 in algorithms:
                start_time = time.time()
                toom2_toom4_toom3(f, g)
                alg2time[243][degree] += time.time() - start_time
            if 52 in algorithms:
                start_time = time.time()
                toom5_toom2(f, g)
                alg2time[52][degree] += time.time() - start_time
            if 53 in algorithms:
                start_time = time.time()
                toom5_toom3(f, g)
                alg2time[53][degree] += time.time() - start_time
            if 54 in algorithms:
                start_time = time.time()
                toom5_toom4(f, g)
                alg2time[54][degree] += time.time() - start_time
            if 522 in algorithms:
                start_time = time.time()
                toom5_toom2_toom2(f, g)
                alg2time[522][degree] += time.time() - start_time
            if 55 in algorithms:
                start_time = time.time()
                toom5_toom5(f, g)
                alg2time[55][degree] += time.time() - start_time
            if 62 in algorithms:
                start_time = time.time()
                toom6_toom2(f, g)
                alg2time[62][degree] += time.time() - start_time
            if 64 in algorithms:
                start_time = time.time()
                toom6_toom4(f, g)
                alg2time[64][degree] += time.time() - start_time
            if 322 in algorithms:
                start_time = time.time()
                toom3_toom2_toom2(f, g)
                alg2time[322][degree] += time.time() - start_time
            if 332 in algorithms:
                start_time = time.time()
                toom3_toom3_toom2(f, g)
                alg2time[332][degree] += time.time() - start_time
            if 333 in algorithms:
                start_time = time.time()
                toom3_toom3_toom3(f, g)
                alg2time[333][degree] += time.time() - start_time
            if 3222 in algorithms:
                start_time = time.time()
                toom3_toom2_toom2_toom2(f, g)
                alg2time[3222][degree] += time.time() - start_time
            if 8 in algorithms:
                start_time = time.time()
                toom8(f, g)
                alg2time[8][degree] += time.time() - start_time
            if 7 in algorithms:
                start_time = time.time()
                toom7(f, g)
                alg2time[7][degree] += time.time() - start_time
            if 72 in algorithms:
                start_time = time.time()
                toom7_toom2(f, g)
                alg2time[72][degree] += time.time() - start_time
            if 73 in algorithms:
                start_time = time.time()
                toom7_toom3(f, g)
                alg2time[73][degree] += time.time() - start_time
            if 74 in algorithms:
                start_time = time.time()
                toom7_toom4(f, g)
                alg2time[74][degree] += time.time() - start_time
            if 82 in algorithms:
                start_time = time.time()
                toom8_toom2(f, g)
                alg2time[82][degree] += time.time() - start_time
            if 83 in algorithms:
                start_time = time.time()
                toom8_toom3(f, g)
                alg2time[83][degree] += time.time() - start_time
            if 65 in algorithms:
                start_time = time.time()
                toom6_toom5(f, g)
                alg2time[65][degree] += time.time() - start_time
            if 532 in algorithms:
                start_time = time.time()
                toom5_toom3_toom2(f, g)
                alg2time[532][degree] += time.time() - start_time
            if 722 in algorithms:
                start_time = time.time()
                toom7_toom2_toom2(f, g)
                alg2time[722][degree] += time.time() - start_time
            if 93 in algorithms:
                start_time = time.time()
                toom9_toom3(f, g)
                alg2time[93][degree] += time.time() - start_time
                
        # progress bar
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
    
    for algorithm in algorithms:
        alg2time[algorithm] = [alg2time[algorithm][degree]/num_trials 
                               for degree in degrees]
    name = {'sb':'Schoolbook', 2:'Toom-2', 3:'Toom-3', 4:'Toom-4',
            5:'Toom-5', 6:'Toom-6', 22:'Toom-2-2', 32:'Toom-3-2',
            432:'Toom-4-3-2', 55:'Toom-5-5', 64:'Toom-6-4',
            3222:'Toom-3-2-2-2', 42:'Toom-4-2', 222:'Toom-2-2-2',
            33: 'Toom-3-3', 8:'Toom-8', 7:'Toom-7', 52:'Toom-5-2',
            43: 'Toom-4-3', 62:'Toom-6-2', 72:'Toom-7-2',
            53:'Toom-5-3', 2222:'Toom-2-2-2-2', 422:'Toom-4-2-2',
            44:'Toom-4-4', 82:'Toom-8-2', 322:'Toom-3-2-2',
            332:'Toom-3-3-2', 54:'Toom-5-4', 522:'Toom-5-2-2',
            73:'Toom-7-3', 83:'Toom-8-3', 243:'Toom-2-4-3',
            74:'Toom-7-4', 333:'Toom-3-3-3', 65:'Toom-6-5',
            532:'Toom-5-3-2', 722:'Toom-7-2-2', 93:'Toom-9-3'}
    for algorithm in algorithms:
        plt.plot(degrees, alg2time[algorithm], lw=3, label=name[algorithm])
    plt.legend(fontsize=13)
    plt.xlabel("Degree", size=15)
    plt.ylabel("Time (seconds)", size=15)
    plt.title("Average Time for Polynomial Multiplication", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    return

def fastest_algorithm(degree, num_trials=100,
                      candidates=('sb', 2, 3, 4, 5, 6, 22, 32, 
                                  222, 42, 33, 8, 7, 52, 43, 62,
                                  72, 2222, 422, 44, 82, 322, 332,
                                  54, 522, 73, 83),
                      progress_bar=False):
    """ This finds which algorithm is fastest for the given degree by taking
        the average over a bunch of trials."""
    alg2time = {}
    for algorithm in candidates:
        alg2time[algorithm] = 0

    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2048, degree)]
        g = [int(x) for x in np.random.randint(0, 2048, degree)]
        if 'sb' in candidates:
            start_time = time.time()
            schoolbook(f, g)
            alg2time['sb'] += time.time() - start_time
        if 2 in candidates:
            start_time = time.time()
            toom2(f, g)
            alg2time[2] += time.time() - start_time
        if 3 in candidates:
            start_time = time.time()
            toom3(f, g)
            alg2time[3] += time.time() - start_time
        if 4 in candidates:
            start_time = time.time()
            toom4(f, g)
            alg2time[4] += time.time() - start_time
        if 5 in candidates:
            start_time = time.time()
            toom5(f, g)
            alg2time[5] += time.time() - start_time
        if 6 in candidates:
            start_time = time.time()
            toom6(f, g)
            alg2time[6] += time.time() - start_time
        if 22 in candidates:
            start_time = time.time()
            toom2_toom2(f, g)
            alg2time[22] += time.time() - start_time
        if 32 in candidates:
            start_time = time.time()
            toom3_toom2(f, g)
            alg2time[32] += time.time() - start_time
        if 33 in candidates:
            start_time = time.time()
            toom3_toom3(f, g)
            alg2time[33] += time.time() - start_time
        if 42 in candidates:
            start_time = time.time()
            toom4_toom2(f, g)
            alg2time[42] += time.time() - start_time
        if 43 in candidates:
            start_time = time.time()
            toom4_toom3(f, g)
            alg2time[43] += time.time() - start_time
        if 52 in candidates:
            start_time = time.time()
            toom5_toom2(f, g)
            alg2time[52] += time.time() - start_time
        if 62 in candidates:
            start_time = time.time()
            toom6_toom2(f, g)
            alg2time[62] += time.time() - start_time
        if 222 in candidates:
            start_time = time.time()
            toom2_toom2_toom2(f, g)
            alg2time[222] += time.time() - start_time
        if 8 in candidates:
            start_time = time.time()
            toom8(f, g)
            alg2time[8] += time.time() - start_time
        if 7 in candidates:
            start_time = time.time()
            toom7(f, g)
            alg2time[7] += time.time() - start_time
        if 72 in candidates:
            start_time = time.time()
            toom7_toom2(f, g)
            alg2time[72] += time.time() - start_time
        if 73 in candidates:
            start_time = time.time()
            toom7_toom3(f, g)
            alg2time[73] += time.time() - start_time
        if 53 in candidates:
            start_time = time.time()
            toom5_toom3(f, g)
            alg2time[53] += time.time() - start_time
        if 2222 in candidates:
            start_time = time.time()
            toom2_toom2_toom2_toom2(f, g)
            alg2time[2222] += time.time() - start_time
        if 322 in candidates:
            start_time = time.time()
            toom3_toom2_toom2(f, g)
            alg2time[322] += time.time() - start_time
        if 422 in candidates:
            start_time = time.time()
            toom4_toom2_toom2(f, g)
            alg2time[422] += time.time() - start_time
        if 44 in candidates:
            start_time = time.time()
            toom4_toom4(f, g)
            alg2time[44] += time.time() - start_time
        if 82 in candidates:
            start_time = time.time()
            toom8_toom2(f, g)
            alg2time[82] += time.time() - start_time
        if 332 in candidates:
            start_time = time.time()
            toom3_toom3_toom2(f, g)
            alg2time[332] += time.time() - start_time
        if 54 in candidates:
            start_time = time.time()
            toom5_toom4(f, g)
            alg2time[54] += time.time() - start_time
        if 522 in candidates:
            start_time = time.time()
            toom5_toom2_toom2(f, g)
            alg2time[522] += time.time() - start_time
        if 83 in candidates:
            start_time = time.time()
            toom8_toom3(f, g)
            alg2time[83] += time.time() - start_time
        # progress bar
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            if progress_bar:
                print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    if progress_bar:
        print("[----------]")
    
    for algorithm in candidates:
        alg2time[algorithm] = alg2time[algorithm]/num_trials
        
    name = {'sb':'SB', 2:'2', 3:'3', 4:'4', 5:'5', 6:'6', 22:'2-2',
            32:'3-2', 42:'4-2', 222:'2-2-2', 33:'3-3', 8:'8', 7:'7',
            52:'5-2', 43:'4-3', 62:'6-2', 72:'7-2', 53:'5-3',
            2222:'2-2-2-2', 422:'4-2-2', 44:'4-4', 82:'8-2',
            322:'3-2-2', 332:'3-2-2', 54:'5-4', 522:'5-2-2',
            73:'7-3', 83:'8-3'}
    return name[min(alg2time, key=alg2time.get)]


def fastest_algorithms(num_trials=500, min_deg=1, max_deg=100,
                       algorithms=('sb', 2, 3, 4, 5, 6, 22, 32,
                                   222, 42, 33, 8, 7, 52, 43, 62,
                                   72, 53, 2222, 422, 44, 82, 322,
                                   332, 54, 522, 73, 83)):
    """ Returns a dictionary mapping each algorithm to a list of the degrees
        where it is the fastest algorithm"""
    name = {'sb':'SB', 2:'2', 3:'3', 4:'4', 5:'5', 6:'6', 22:'2-2',
            32:'3-2', 42:'4-2', 222:'2-2-2', 33:'3-3', 8:'8', 7:'7',
            52:'5-2', 43:'4-3', 62:'6-2', 72:'7-2', 53:'5-3',
            2222:'2-2-2-2', 422:'4-2-2', 44:'4-4', 82:'8-2',
            322:'3-2-2', 332:'3-3-2', 54:'5-4', 522:'5-2-2',
            73:'7-3', 83:'8-3'}
    alg2degrees = {name[algorithm]:[] for algorithm in algorithms}
    for degree in range(min_deg, max_deg+1):
        winner = fastest_algorithm(degree, num_trials=num_trials,
                                   candidates=algorithms)
        alg2degrees[winner].append(degree)
        print(degree)
    for algorithm in alg2degrees:
        int_list = alg2degrees[algorithm]
        print(algorithm.ljust(5) + ": " + list_to_intervals(int_list))
    return alg2degrees

def list_to_intervals(L):
    """ This converts an increasing list of distinct positive integers into
        a string that uses a dash for a run of consecutives. 
        For example, [1,2,3,5,7,8] becomes "1-3, 5, 7-8" """
    if len(L) == 0:
        return ""
    if len(L) == 1:
        return str(L[0])
    
    s = str(L[0])
    prev_int = L[0]
    cur_run_len = 0
    for i in L[1:]:
        if cur_run_len == 0 and i == prev_int + 1:
            s += '-'
            prev_int = i
            cur_run_len = 1
        elif cur_run_len > 0 and i == prev_int + 1:
            prev_int = i
            cur_run_len += 1
        elif cur_run_len > 0 and i != prev_int + 1:
            s += str(prev_int)
            s += ", "
            s += str(i)
            prev_int = i
            cur_run_len = 0
        elif cur_run_len == 0 and i != prev_int + 1:
            s += ", "
            s += str(i)
            prev_int = i
    if s[-1] == '-':     # a cheap fix
        s += str(L[-1])
    return s

    
#fastest_algorithms(num_trials=1500, min_deg=131, max_deg=131,
#                   algorithms=(5, 6, 72, 53, 2222, 422, 44, 82,
#                                   33, 8, 7, 52, 43, 62, 322))




 
#fastest_algorithms(num_trials=50, min_deg=819, max_deg=822,
 #                  algorithms=(333, 74, 55)) 

plot_algorithms([333,74,532,65,93,722], num_trials=2000, min_deg=816, max_deg=825,
                filename="compare_6_at_821.jpg")
#plot_algorithms([55, 64, 83], num_trials=1000, min_deg=735, max_deg=750, filename="best_743.jpg")
#plot_algorithms([432, 55, 64, 'sb'], num_trials=300, min_deg=735, max_deg=750,
#                filename='degree742compare4.jpg')  

#plot_algorithms([42, 6], num_trials=500, min_deg=120, max_deg=300,
#                filename='toom42toom6threshold.jpg')   
#plot_algorithms([42, 32], num_trials=500, min_deg=20, max_deg=200,
#plot_algorithms([42, 222], num_trials=500, min_deg=20, max_deg=200,
#               filename='toom222toom42threshold.jpg')  
#plot_algorithms([32, 222], num_trials=500, min_deg=80, max_deg=300,
#                filename='toom222toom32threshold.jpg')       
                
      
                
                