# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:48:10 2019

@author: Marcus
"""

import numpy as np
import matplotlib.pyplot as plt
import time


def format_time(x):
    if x < 60:
        return str(int(x)) + " seconds"
    elif x < 3600:
        return str(int(x//60)) + "." + str(int((x-60*(x//60))/60*10)) + " minutes"
    else:
        return str(int(x//3600)) + "." + str(int((x-3600*(x//3600))/3600*10)) + " hours"



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

def toom4_schoolbook(f, g):
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

#    r0 = r[0]
#    r6 = r['infinity']
#    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r6[i]
#          for i in range(len(r0))]
#    r4 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 64*r6[i])//4 - T1[i])//3
#          for i in range(len(r0))]
#    r2 = [T1[i] - r4[i] for i in range(len(r0))]
#    T2 = [(r[1][i] - r[-1][i])//2 for i in range(len(r0))]
#    T3 = [((r[2][i] - r[-2][i])//4 - T2[i])//3 for i in range(len(r0))]
#    T4 = [((r[3][i] - r0[i] - 9*r2[i] - 81*r4[i] - 729*r6[i])//3 - T2[i])//8
#          for i in range(len(r0))]
#    r3 = [2*T3[i] - T4[i] for i in range(len(r0))]
#    r5 = [(T3[i] - r3[i])//10 for i in range(len(r0))]
#    r1 = [T2[i] - r3[i] - r5[i] for i in range(len(r0))]
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


def optimized_toom_6(f, g):
    """ Do the 6-4 algorithm"""
    n = 6
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    r = {eval_list[i] : toom4_schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
    r10 = r['infinity']
    T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i]
          for i in range(len(r0))]
    T2 = [(((r[2][i] + r[-2][i])//2 - r0[i] - 1024*r10[i])//4 - T1[i])//3
          for i in range(len(r0))]
    T3 = [(((r[3][i] + r[-3][i])//2 - r0[i] - 59049*r10[i])//9 - T1[i])//8
          for i in range(len(r0))]
    T4 = [(((r[4][i] + r[-4][i])//2 - r0[i] - 1048576*r10[i])//16 - T1[i])//15
          for i in range(len(r0))]
    T5 = [(T3[i] - T2[i])//5 for i in range(len(r0))]
    T6 = [(T4[i] - T3[i])//7 for i in range(len(r0))]
    r8 = [(T6[i] - T5[i])//12 for i in range(len(r0))]
    r6 = [T5[i] - 14*r8[i] for i in range(len(r0))]
    r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in range(len(r0))]
    r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in range(len(r0))]
    
    T7 = [(r[1][i] - r[-1][i])//2 for i in range(len(r0))]
    T8 = [((r[2][i] - r[-2][i])//4 - T7[i])//3 for i in range(len(r0))]
    T9 = [((r[3][i] - r[-3][i])//6 - T7[i])//8 for i in range(len(r0))]
    T10 = [((r[4][i] - r[-4][i])//8 - T7[i])//15 for i in range(len(r0))]
    T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] - 15625*r6[i] 
    - 390625*r8[i] - 9765625*r10[i])//5 - T7[i])//24 for i in range(len(r0))]
    T12 = [(T9[i] - T8[i])//5 for i in range(len(r0))]
    T13 = [(T10[i] - T9[i])//7 for i in range(len(r0))]
    T14 = [(T11[i] - T10[i])//9 for i in range(len(r0))]
    T15 = [(T13[i] - T12[i])//12 for i in range(len(r0))]
    T16 = [(T14[i] - T13[i])//16 for i in range(len(r0))]
    r9 = [(T16[i] - T15[i])//21 for i in range(len(r0))]
    r7 = [T15[i] - 30*r9[i] for i in range(len(r0))]
    r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in range(len(r0))]
    r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in range(len(r0))]
    r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in range(len(r0))]

    for i in range(len(r0)):
        print(i, r1[i] % 2048)

    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]


def test_optimized_toom_6():
    for poly_length in range(100, 110):
        for _ in range(50):
            f = [int(x) for x in np.random.randint(0, 2048, size=poly_length)]
            g = [int(x) for x in np.random.randint(0, 2048, size=poly_length)]
            exp = schoolbook(f, g)
            obs = optimized_toom_6(f, g)
            for i in range(len(exp)):
                if exp[i] != obs[i]:
                    print(f, g)
                    for j in range(len(exp)):
                        print(j, exp[j], obs[j], exp[j] - obs[j])
                    return

#test_optimized_toom_6()
f = range(1, 769)
g = [768 - i for i in range(768)]
optimized_toom_6(f,g)