# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:25:59 2019

@author: Marcus
"""
import numpy as np
from scipy import stats
# =========================================
#
#           Helper Functions
#
# ========================================

def schoolbook_mod(f, g, m):
    """ Uses schoolbook multiplication to multiply f and g mod 
        m. Returns the product as a list"""
    d = len(f) + len(g) - 1
    
    # initialize a list of zeros
    product = [0]*d
    
    # distribute through all possible combinations of coefficients
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] = (product[i+j] + f[i]*g[j]) % m
    return product


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

def evaluate_blocks_mod(blocks, value, m):
    """ blocks is a list of lists, each list is the coefficients of a
        polynomial. But each list a coefficient. For example, if blocks is
        [[1,2],[3,4],[5,6]] and value is -2, we return
        [1,2] + [-6,-8] + [20,24] = [15, 18].  If the value is infinity,
        we return the leading coefficient.
        This does it mod m"""
        
    if value == 'infinity':
        return blocks[-1]
    
    # initialize an empty list of the right length
    answer = [0]*len(blocks[0])
    
    coefficient = 1
    for i in range(len(blocks)):
        for j in range(len(blocks[0])):
            answer[j] = (answer[j] + coefficient*blocks[i][j]) % m
        coefficient = (coefficient*value) % m
    return answer

def evaluate_blocks_list_mod(blocks, values, m):
    """ Evaluates the blocks on a list of values, and returns a list"""
    answer = []
    for value in values:
        answer.append(evaluate_blocks_mod(blocks, value, m))
    return answer


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



def max_power_of_two(a, b, max_n=32):
    """ Returns the largest power of two such that a = b mod 2^n
        up to 2**max_n"""
    n = max_n
    while (a % (2**n)) != (b % (2**n)):
        n -= 1
    return n

def max_power_of_two_list(A, B, max_n):
    min_so_far = max_n
    for i in range(len(A)):
        if max_power_of_two(A[i], B[i], max_n=max_n) < min_so_far:
            min_so_far = max_power_of_two(A[i], B[i], max_n=max_n)
            #print(A[i], B[i], min_so_far)
    return min_so_far

# ========================================
#
#          The New Function
#
# =======================================

def check_toom4_cost(num_trials=50):
    exp = 16
    m = 2**exp
    inv3 = 43691
    inv15 = 61167
    eval_list = (0, 1, -1, 2, -2, 3, 'infinity')
    # this dictionary maps each variable we care about to a list of
    # the costs at that variable from each trial
    D = {}
    for i in eval_list:
        label_string = "r({})".format(i)
        D[label_string] = []
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        D[label_string] = []
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2048, size=16)]
        g = [int(x) for x in np.random.randint(0, 2048, size=16)]
        fblocks = split(f, 4)
        gblocks = split(g, 4)
        f_eval = evaluate_blocks_list(fblocks, eval_list)
        g_eval = evaluate_blocks_list(gblocks, eval_list)
        f_eval_m = evaluate_blocks_list_mod(fblocks, eval_list, m)
        g_eval_m = evaluate_blocks_list_mod(gblocks, eval_list, m)
        # perform the recursive multiplication
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
             for i in range(len(f_eval))}
        r_m = {eval_list[i] : schoolbook_mod(f_eval_m[i], g_eval_m[i], m)
             for i in range(len(f_eval_m))}
        
        r0 = r[0]
        r0_m = r_m[0]
        
        r6 = r['infinity']
        r6_m = r_m['infinity']
        
        r4 = [(6*r0[i] 
              - 120*r6[i] 
              + (r[2][i]+r[-2][i]) 
              - 4*(r[1][i]+r[-1][i])
              )//24 
            for i in range(len(r0))]
        r4_m = [6*r0_m[i] % m for i in range(len(r0_m))]
        for i in range(len(r0)):
            r4_m[i] = (r4_m[i] - 120*r6_m[i]) % m
            r4_m[i] = (r4_m[i] + r_m[2][i] + r_m[-2][i]) % m
            r4_m[i] = (r4_m[i] - 4*(r_m[1][i] + r_m[-1][i])) % m
            r4_m[i] = r4_m[i] // 8
            r4_m[i] = (r4_m[i] * inv3) % m
        
        r2 = [(r[1][i] 
              + r[-1][i])//2 
              - r0[i] 
              - r4[i] 
              - r6[i] 
              for i in range(len(r0))]
        r2_m = [r_m[1][i] + r_m[-1][i] for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r2_m[i] = (r2_m[i] - 2*r0_m[i]) % m
            r2_m[i] = (r2_m[i] - 2*r4_m[i]) % m
            r2_m[i] = (r2_m[i] - 2*r6_m[i]) % m
            r2_m[i] = r2_m[i] // 2
        
        r5 = [(r[3][i] 
              - 4*r[2][i] 
              + 5*r[1][i] 
              - 2*r0[i] 
              + 2*r2[i] 
              - 22*r4[i] 
              - 478*r6[i])//120 
              for i in range(len(r0))]
        r5_m = [r_m[3][i] for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r5_m[i] = (r5_m[i] - 4*r_m[2][i]) % m
            r5_m[i] = (r5_m[i] + 5*r_m[1][i]) % m
            r5_m[i] = (r5_m[i] - 2*r0_m[i]) % m
            r5_m[i] = (r5_m[i] + 2*r2_m[i]) % m
            r5_m[i] = (r5_m[i] - 22*r4_m[i]) % m
            r5_m[i] = (r5_m[i] - 478*r6_m[i]) % m
            r5_m[i] = r5_m[i] // 8
            r5_m[i] = (r5_m[i] * inv15) % m
        
        r3 = [(r[2][i] 
              - 2*r[1][i] 
              + r0[i] 
              - 2*r2[i] 
              - 14*r4[i] 
              - 30*r5[i] 
              - 62*r6[i])//6 
             for i in range(len(r0))]
        r3_m = [r_m[2][i] for i in range(len(r0_m))]
        for i in range(len(r0)):
            r3_m[i] = (r3_m[i] - 2*r_m[1][i]) % m
            r3_m[i] = (r3_m[i] + r0_m[i]) % m
            r3_m[i] = (r3_m[i] - 2*r2_m[i]) % m
            r3_m[i] = (r3_m[i] - 14*r4_m[i]) % m
            r3_m[i] = (r3_m[i] - 30*r5_m[i]) % m
            r3_m[i] = (r3_m[i] - 62*r6_m[i]) % m
            r3_m[i] = r3_m[i] // 2
            r3_m[i] = (r3_m[i] * inv3) % m
        
        r1 = [r[1][i] 
              - r0[i] 
              - r2[i] 
              - r3[i] 
              - r4[i] 
              - r5[i] 
              - r6[i] 
              for i in range(len(r0))]
        r1_m = [r_m[1][i] - r0_m[i] for i in range(len(r0_m))]
        for i in range(len(r0)):
            r1_m[i] = (r1_m[i] - r2_m[i]) % m
            r1_m[i] = (r1_m[i] - r3_m[i]) % m
            r1_m[i] = (r1_m[i] - r4_m[i]) % m
            r1_m[i] = (r1_m[i] - r5_m[i]) % m
            r1_m[i] = (r1_m[i] - r6_m[i]) % m 
        
        r_coefs = (r0, r1, r2, r3, r4, r5, r6)
        r_coefs_m = (r0_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m)
        for i in eval_list:
            label_string = "r({})".format(i)
            D[label_string].append(exp - max_power_of_two_list(r[i], r_m[i], exp))
        for i in range(2*4 - 1):
            label_string = "r{}".format(i)
            D[label_string].append(exp - max_power_of_two_list(r_coefs[i], r_coefs_m[i], exp))
    for i in eval_list:
        label_string = "r({})".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    return

def check_optimized_toom4_cost(num_trials=50):
    exp = 16
    m = 2**exp
    inv3 = 43691
    inv5 = 52429
    eval_list = (0, 1, -1, 2, -2, 3, 'infinity')
    # this dictionary maps each variable we care about to a list of
    # the costs at that variable from each trial
    D = {}
    for i in eval_list:
        label_string = "r({})".format(i)
        D[label_string] = []
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        D[label_string] = []
    for i in range(1, 4+1):
        label_string = "T{}".format(i)
        D[label_string] = []
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2048, size=16)]
        g = [int(x) for x in np.random.randint(0, 2048, size=16)]
        fblocks = split(f, 4)
        gblocks = split(g, 4)
        f_eval = evaluate_blocks_list(fblocks, eval_list)
        g_eval = evaluate_blocks_list(gblocks, eval_list)
        f_eval_m = evaluate_blocks_list_mod(fblocks, eval_list, m)
        g_eval_m = evaluate_blocks_list_mod(gblocks, eval_list, m)
        # perform the recursive multiplication
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
             for i in range(len(f_eval))}
        r_m = {eval_list[i] : schoolbook_mod(f_eval_m[i], g_eval_m[i], m)
             for i in range(len(f_eval_m))}
        
        r0 = r[0]
        r0_m = r_m[0]
        
        r6 = r['infinity']
        r6_m = r_m['infinity']
       
        T1 = [(r[1][i] + r[-1][i]) // 2 
              - r0[i] 
              - r6[i] for i in range(len(r0))]
        T1_m = [(r_m[1][i] + r_m[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T1_m[i] = T1_m[i] // 2
            T1_m[i] = (T1_m[i] - r0_m[i]) % m
            T1_m[i] = (T1_m[i] - r6_m[i]) % m
            
        r4 = [(((r[2][i] + r[-2][i]) // 2
              - r0[i]
              - 64*r6[i]) // 4
              - T1[i]) // 3
            for i in range(len(r0))]
        r4_m = [(r_m[2][i] + r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r4_m[i] = r4_m[i] // 2
            r4_m[i] = (r4_m[i] - r0_m[i]) % m
            r4_m[i] = (r4_m[i] - 64*r6_m[i]) % m
            r4_m[i] = r4_m[i] // 4
            r4_m[i] = (r4_m[i] - T1_m[i]) % m
            r4_m[i] = (r4_m[i] * inv3) % m
        
        r2 = [T1[i] - r4[i] for i in range(len(r0))]
        r2_m = [(T1_m[i] - r4_m[i]) % m for i in range(len(r0_m))]
        
        T2 = [(r[1][i] - r[-1][i]) // 2 for i in range(len(r0))]
        T2_m = [(r_m[1][i] - r[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T2_m[i] = T2_m[i] // 2
        
        T3 = [((r[2][i] - r[-2][i]) // 4 - T2[i]) //3 for i in range(len(r0))]
        T3_m = [(r_m[2][i] - r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T3_m[i] = T3_m[i] // 4
            T3_m[i] = (T3_m[i] - T2_m[i]) % m
            T3_m[i] = (T3_m[i] * inv3) % m
        
        T4 = [((r[3][i]
               - r0[i]
               - 9*r2[i]
               - 81*r4[i]
               - 729*r6[i]) // 3
               - T2[i]) // 8
            for i in range(len(r0))]
        T4_m = [(r_m[3][i] - r0_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T4_m[i] = (T4_m[i] - 9*r2_m[i]) % m
            T4_m[i] = (T4_m[i] - 81*r4_m[i]) % m
            T4_m[i] = (T4_m[i] - 729*r6_m[i]) % m
            T4_m[i] = (T4_m[i] * inv3) % m
            T4_m[i] = (T4_m[i] - T2_m[i]) % m
            T4_m[i] = T4_m[i] // 8
        
        r5 = [(T4[i] - T3[i]) // 5 for i in range(len(r0))]
        r5_m = [(T4_m[i] - T3_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r5_m[i] = (r5_m[i] * inv5) % m
            
        r3 = [T3[i] - 5*r5[i] for i in range(len(r0))]
        r3_m = [(T3_m[i] - 5*r5_m[i]) % m for i in range(len(r0_m))]
        
        r1 = [T2[i] - r3[i] - r5[i] for i in range(len(r0))]
        r1_m = [(T2_m[i] - r3_m[i] - r5_m[i]) % m for i in range(len(r0_m))]
            
        r_coefs = (r0, r1, r2, r3, r4, r5, r6)
        r_coefs_m = (r0_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m)
        
        T_coefs = (T1, T2, T3, T4)
        T_coefs_m = (T1_m, T2_m, T3_m, T4_m)
        
        for i in eval_list:
            label_string = "r({})".format(i)
            D[label_string].append(exp - max_power_of_two_list(r[i], r_m[i], exp))
        for i in range(2*4 - 1):
            label_string = "r{}".format(i)
            D[label_string].append(exp - max_power_of_two_list(r_coefs[i], r_coefs_m[i], exp))
        for i in range(1, 5):
             label_string = "T{}".format(i)
             D[label_string].append(exp - max_power_of_two_list(T_coefs[i-1], T_coefs_m[i-1], exp))
    for i in eval_list:
        label_string = "r({})".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(1,5):
        label_string = "T{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    return      
            
def check_other_toom4_cost(num_trials=50):
    exp = 16
    m = 2**exp
    inv3 = 43691
    inv5 = 52429
    eval_list = (0, 1, -1, 2, -2, 3, 'infinity')
    # this dictionary maps each variable we care about to a list of
    # the costs at that variable from each trial
    D = {}
    for i in eval_list:
        label_string = "r({})".format(i)
        D[label_string] = []
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        D[label_string] = []
    for i in range(1, 4+1):
        label_string = "T{}".format(i)
        D[label_string] = []
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2048, size=16)]
        g = [int(x) for x in np.random.randint(0, 2048, size=16)]
        fblocks = split(f, 4)
        gblocks = split(g, 4)
        f_eval = evaluate_blocks_list(fblocks, eval_list)
        g_eval = evaluate_blocks_list(gblocks, eval_list)
        f_eval_m = evaluate_blocks_list_mod(fblocks, eval_list, m)
        g_eval_m = evaluate_blocks_list_mod(gblocks, eval_list, m)
        # perform the recursive multiplication
        C = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
             for i in range(len(f_eval))}
        C_m = {eval_list[i] : schoolbook_mod(f_eval_m[i], g_eval_m[i], m)
             for i in range(len(f_eval_m))}
        
        C0 = C[0]
        C0_m = C_m[0]
        
        C6 = C['infinity']
        C6_m = C_m['infinity']
       
        V0 = [(C[1][i] + C[-1][i]) // 2 
              - C0[i] 
              - C6[i] for i in range(len(C0))]
        V0_m = [(C_m[1][i] + C_m[-1][i]) % m for i in range(len(C0_m))]
        for i in range(len(C0_m)):
            V0_m[i] = V0_m[i] // 2
            V0_m[i] = (V0_m[i] - C0_m[i]) % m
            V0_m[i] = (V0_m[i] - C6_m[i]) % m
            
        V1 = [(C[2][i] + C[-2][i] - 2*C0[i])]
            
        r4 = [(((r[2][i] + r[-2][i]) // 2
              - r0[i]
              - 64*r6[i]) // 4
              - T1[i]) // 3
            for i in range(len(r0))]
        r4_m = [(r_m[2][i] + r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r4_m[i] = r4_m[i] // 2
            r4_m[i] = (r4_m[i] - r0_m[i]) % m
            r4_m[i] = (r4_m[i] - 64*r6_m[i]) % m
            r4_m[i] = r4_m[i] // 4
            r4_m[i] = (r4_m[i] - T1_m[i]) % m
            r4_m[i] = (r4_m[i] * inv3) % m
        
        r2 = [T1[i] - r4[i] for i in range(len(r0))]
        r2_m = [(T1_m[i] - r4_m[i]) % m for i in range(len(r0_m))]
        
        T2 = [(r[1][i] - r[-1][i]) // 2 for i in range(len(r0))]
        T2_m = [(r_m[1][i] - r[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T2_m[i] = T2_m[i] // 2
        
        T3 = [((r[2][i] - r[-2][i]) // 4 - T2[i]) //3 for i in range(len(r0))]
        T3_m = [(r_m[2][i] - r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T3_m[i] = T3_m[i] // 4
            T3_m[i] = (T3_m[i] - T2_m[i]) % m
            T3_m[i] = (T3_m[i] * inv3) % m
        
        T4 = [((r[3][i]
               - r0[i]
               - 9*r2[i]
               - 81*r4[i]
               - 729*r6[i]) // 3
               - T2[i]) // 8
            for i in range(len(r0))]
        T4_m = [(r_m[3][i] - r0_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T4_m[i] = (T4_m[i] - 9*r2_m[i]) % m
            T4_m[i] = (T4_m[i] - 81*r4_m[i]) % m
            T4_m[i] = (T4_m[i] - 729*r6_m[i]) % m
            T4_m[i] = (T4_m[i] * inv3) % m
            T4_m[i] = (T4_m[i] - T2_m[i]) % m
            T4_m[i] = T4_m[i] // 8
        
        r5 = [(T4[i] - T3[i]) // 5 for i in range(len(r0))]
        r5_m = [(T4_m[i] - T3_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r5_m[i] = (r5_m[i] * inv5) % m
            
        r3 = [T3[i] - 5*r5[i] for i in range(len(r0))]
        r3_m = [(T3_m[i] - 5*r5_m[i]) % m for i in range(len(r0_m))]
        
        r1 = [T2[i] - r3[i] - r5[i] for i in range(len(r0))]
        r1_m = [(T2_m[i] - r3_m[i] - r5_m[i]) % m for i in range(len(r0_m))]
            
        r_coefs = (r0, r1, r2, r3, r4, r5, r6)
        r_coefs_m = (r0_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m)
        
        T_coefs = (T1, T2, T3, T4)
        T_coefs_m = (T1_m, T2_m, T3_m, T4_m)
        
        for i in eval_list:
            label_string = "r({})".format(i)
            D[label_string].append(exp - max_power_of_two_list(r[i], r_m[i], exp))
        for i in range(2*4 - 1):
            label_string = "r{}".format(i)
            D[label_string].append(exp - max_power_of_two_list(r_coefs[i], r_coefs_m[i], exp))
        for i in range(1, 5):
             label_string = "T{}".format(i)
             D[label_string].append(exp - max_power_of_two_list(T_coefs[i-1], T_coefs_m[i-1], exp))
    for i in eval_list:
        label_string = "r({})".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(2*4 - 1):
        label_string = "r{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(1,5):
        label_string = "T{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    return

def check_optimized_toom5_cost(num_trials=50):
    exp = 32
    m = 2**exp
    inv3 = 2863311531
    inv5 = 3435973837
    inv7 = 3067833783
    inv9 = 954437177
    inv15 = 4008636143
    eval_list = (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    # this dictionary maps each variable we care about to a list of
    # the costs at that variable from each trial
    D = {}
    for i in eval_list:
        label_string = "r({})".format(i)
        D[label_string] = []
    for i in range(2*5 - 1):
        label_string = "r{}".format(i)
        D[label_string] = []
    for i in range(1, 9+1):
        label_string = "T{}".format(i)
        D[label_string] = []
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2**(20), size=25)]
        g = [int(x) for x in np.random.randint(0, 2**(20), size=25)]
        fblocks = split(f, 5)
        gblocks = split(g, 5)
        f_eval = evaluate_blocks_list(fblocks, eval_list)
        g_eval = evaluate_blocks_list(gblocks, eval_list)
        f_eval_m = evaluate_blocks_list_mod(fblocks, eval_list, m)
        g_eval_m = evaluate_blocks_list_mod(gblocks, eval_list, m)
        # perform the recursive multiplication
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
             for i in range(len(f_eval))}
        r_m = {eval_list[i] : schoolbook_mod(f_eval_m[i], g_eval_m[i], m)
             for i in range(len(f_eval_m))}
        
        r0 = r[0]
        r0_m = r_m[0]

        r8 = r['infinity']
        r8_m = r_m['infinity']

        T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r8[i] for i in range(len(r0))]
        T1_m = [(r_m[1][i] + r_m[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T1_m[i] = T1_m[i] // 2
            T1_m[i] = (T1_m[i] - r0_m[i]) % m
            T1_m[i] = (T1_m[i] - r8_m[i]) % m
            
    
        T2 = [(((r[2][i] + r[-2][i])//2 
              - r0[i] 
              - 256*r8[i])//4
            - T1[i]) // 3 for i in range(len(r0))]
        T2_m = [(r_m[2][i] + r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T2_m[i] = T2_m[i] // 2
            T2_m[i] = (T2_m[i] - r0_m[i]) % m
            T2_m[i] = (T2_m[i] - 256*r8_m[i]) % m
            T2_m[i] = T2_m[i] // 4
            T2_m[i] = (T2_m[i] - T1_m[i]) % m
            T2_m[i] = (T2_m[i] * inv3) % m
    
        T3 = [(((r[3][i] + r[-3][i])//2 
               - r0[i] 
               - 6561*r8[i]) // 9 
             - T1[i]) // 8 for i in range(len(r0))]
        T3_m = [(r_m[3][i] + r_m[-3][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T3_m[i] = T3_m[i] // 2
            T3_m[i] = (T3_m[i] - r0_m[i]) % m
            T3_m[i] = (T3_m[i] - 6561*r8_m[i]) % m
            T3_m[i] = (T3_m[i] * inv9) % m
            T3_m[i] = (T3_m[i] - T1_m[i]) % m
            T3_m[i] = T3_m[i] // 8
    
        r6 = [(T3[i] - T2[i]) // 5 for i in range(len(r0))]
        r6_m = [((T3_m[i] - T2_m[i]) * inv5 ) % m for i in range(len(r0_m))]
    
        r4 = [T2[i] - 5*r6[i] for i in range(len(r0))]
        r4_m = [(T2_m[i] - 5*r6_m[i]) % m for i in range(len(r0_m))]
    
        r2 = [T1[i] - r4[i] - r6[i] for i in range(len(r0))]
        r2_m = [(T1_m[i] - r4_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r2_m[i] = (r2_m[i] - r6_m[i]) % m
    
        T4 = [(r[1][i] - r[-1][i])//2 for i in range(len(r0))]
        T4_m = [(r_m[1][i] - r_m[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T4_m[i] = T4_m[i] // 2
    
        T5 = [((r[2][i] - r[-2][i])//4
              - T4[i]) // 3 for i in range(len(r0))]
        T5_m = [(r_m[2][i] - r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T5_m[i] = T5_m[i] // 4
            T5_m[i] = (T5_m[i] - T4_m[i]) % m
            T5_m[i] = (T5_m[i] * inv3) % m
    
        T6 = [((r[3][i] - r[-3][i])//6
              - T4[i]) // 8 for i in range(len(r0))]
        T6_m = [(r_m[3][i] - r_m[-3][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T6_m[i] = T6_m[i] // 2
            T6_m[i] = (T6_m[i] * inv3) % m
            T6_m[i] = (T6_m[i] - T4_m[i]) % m
            T6_m[i] = T6_m[i] // 8
    
        T7 = [((r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] 
                       - 4096*r6[i] - 65536*r8[i])//4
             - T4[i]) // 15 for i in range(len(r0))]
        T7_m = [(r_m[4][i] - r0_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T7_m[i] = (T7_m[i] - 16*r2_m[i]) % m
            T7_m[i] = (T7_m[i] - 256*r4_m[i]) % m
            T7_m[i] = (T7_m[i] - 4096*r6_m[i]) % m
            T7_m[i] = (T7_m[i] - 65536*r8_m[i]) % m
            T7_m[i] = T7_m[i] // 4
            T7_m[i] = (T7_m[i] - T4_m[i]) % m
            T7_m[i] = (T7_m[i] * inv15) % m
    
        T8 = [(T6[i] - T5[i]) // 5 for i in range(len(r0))]
        T8_m = [((T6_m[i] - T5_m[i]) * inv5) % m for i in range(len(r0_m))]
    
        T9 = [(T7[i] - T6[i]) // 7 for i in range(len(r0))]
        T9_m = [((T7_m[i] - T6_m[i]) * inv7) % m for i in range(len(r0_m))]
               
        r7 = [(T9[i] - T8[i]) // 12 for i in range(len(r0))]
        r7_m = [(T9_m[i] - T8_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r7_m[i] = r7_m[i] // 4
            r7_m[i] = (r7_m[i] * inv3) % m
    
        r5 = [T8[i] - 14*r7[i] for i in range(len(r0))]
        r5_m = [(T8_m[i] - 14*r7_m[i]) % m for i in range(len(r0_m))]
    
        r3 = [T5[i] - 5*r5[i] - 21*r7[i] for i in range(len(r0))]
        r3_m = [(T5_m[i] - 5*r5_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r3_m[i] = (r3_m[i] - 21*r7_m[i]) % m
    
        r1 = [T4[i] - r3[i] - r5[i] - r7[i] for i in range(len(r0))]
        r1_m = [(T4_m[i] - r3_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r1_m[i] = (r1_m[i] - r5_m[i]) % m
            r1_m[i] = (r1_m[i] - r7_m[i]) % m
    
        r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
        r_coefs_m = (r0_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m, r7_m, r8_m)
        
        T_coefs = (T1, T2, T3, T4, T5, T6, T7, T8, T9)
        T_coefs_m = (T1_m, T2_m, T3_m, T4_m, T5_m, T6_m, T7_m, T8_m, T9_m)
        
        for i in eval_list:
            label_string = "r({})".format(i)
            D[label_string].append(exp - max_power_of_two_list(r[i], r_m[i], exp))
        for i in range(2*5 - 1):
            label_string = "r{}".format(i)
            D[label_string].append(exp - max_power_of_two_list(r_coefs[i], r_coefs_m[i], exp))
        for i in range(1, 10):
             label_string = "T{}".format(i)
             D[label_string].append(exp - max_power_of_two_list(T_coefs[i-1], T_coefs_m[i-1], exp))
    for i in eval_list:
        label_string = "r({})".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(2*5 - 1):
        label_string = "r{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(1,10):
        label_string = "T{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    return   


def check_optimized_toom6_cost(num_trials=50):
    exp = 32
    m = 2**exp
    inv3 = 2863311531
    inv5 = 3435973837
    inv7 = 3067833783
    inv9 = 954437177
    inv15 = 4008636143
    inv21 = 1022611261
    eval_list = (0, 1, -1, 2, -2, 3, -3, 4,-4, 5, 'infinity')
    # this dictionary maps each variable we care about to a list of
    # the costs at that variable from each trial
    D = {}
    for i in eval_list:
        label_string = "r({})".format(i)
        D[label_string] = []
    for i in range(2*6 - 1):
        label_string = "r{}".format(i)
        D[label_string] = []
    for i in range(1, 16+1):
        label_string = "T{}".format(i)
        D[label_string] = []
    for _ in range(num_trials):
        f = [int(x) for x in np.random.randint(0, 2**(20), size=25)]
        g = [int(x) for x in np.random.randint(0, 2**(20), size=25)]
        fblocks = split(f, 5)
        gblocks = split(g, 5)
        f_eval = evaluate_blocks_list(fblocks, eval_list)
        g_eval = evaluate_blocks_list(gblocks, eval_list)
        f_eval_m = evaluate_blocks_list_mod(fblocks, eval_list, m)
        g_eval_m = evaluate_blocks_list_mod(gblocks, eval_list, m)
        # perform the recursive multiplication
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
             for i in range(len(f_eval))}
        r_m = {eval_list[i] : schoolbook_mod(f_eval_m[i], g_eval_m[i], m)
             for i in range(len(f_eval_m))}
        
        r0 = r[0]
        r0_m = r_m[0]

        r10 = r['infinity']
        r10_m = r_m['infinity']

        T1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r10[i] for i in range(len(r0))]
        T1_m = [(r_m[1][i] + r_m[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T1_m[i] = T1_m[i] // 2
            T1_m[i] = (T1_m[i] - r0_m[i]) % m
            T1_m[i] = (T1_m[i] - r10_m[i]) % m
            
    
        T2 = [(((r[2][i] + r[-2][i])//2 
              - r0[i] 
              - 1024*r10[i])//4
            - T1[i]) // 3 for i in range(len(r0))]
        T2_m = [(r_m[2][i] + r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T2_m[i] = T2_m[i] // 2
            T2_m[i] = (T2_m[i] - r0_m[i]) % m
            T2_m[i] = (T2_m[i] - 1024*r10_m[i]) % m
            T2_m[i] = T2_m[i] // 4
            T2_m[i] = (T2_m[i] - T1_m[i]) % m
            T2_m[i] = (T2_m[i] * inv3) % m
    
        T3 = [(((r[3][i] + r[-3][i])//2 
               - r0[i] 
               - 59049*r10[i]) // 9 
             - T1[i]) // 8 for i in range(len(r0))]
        T3_m = [(r_m[3][i] + r_m[-3][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T3_m[i] = T3_m[i] // 2
            T3_m[i] = (T3_m[i] - r0_m[i]) % m
            T3_m[i] = (T3_m[i] - 59049*r10_m[i]) % m
            T3_m[i] = (T3_m[i] * inv9) % m
            T3_m[i] = (T3_m[i] - T1_m[i]) % m
            T3_m[i] = T3_m[i] // 8
            
        T4 = [(((r[4][i] + r[-4][i])//2 
               - r0[i] 
               - 1048576*r10[i]) // 16 
             - T1[i]) // 15 for i in range(len(r0))]
        T4_m = [(r_m[4][i] + r_m[-4][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T4_m[i] = T4_m[i] // 2
            T4_m[i] = (T4_m[i] - r0_m[i]) % m
            T4_m[i] = (T4_m[i] - 1048576*r10_m[i]) % m
            T4_m[i] = T4_m[i] // 16
            T4_m[i] = (T4_m[i] - T1_m[i]) % m
            T4_m[i] = (T4_m[i] * inv15) % m
            
        T5 = [(T3[i] - T2[i]) // 5 for i in range(len(r0))]
        T5_m = [((T3_m[i] - T2_m[i]) * inv5) % m for i in range(len(r0_m))]
        
        T6 = [(T4[i] - T3[i]) // 7 for i in range(len(r0))]
        T6_m = [((T4_m[i] - T3_m[i]) * inv7) % m for i in range(len(r0_m))]
        
        r8 = [(T6[i] - T5[i]) // 12 for i in range(len(r0))]
        r8_m = [(T6_m[i] - T5_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r8_m[i] = r8_m[i] // 4
            r8_m[i] = (r8_m[i] * inv3) % m
            
        r6 = [T5[i] - 14*r8[i] for i in range(len(r0))]
        r6_m = [(T5_m[i] - 14*r8_m[i]) % m for i in range(len(r0_m))]
    
        r4 = [T2[i] - 5*r6[i] - 21*r8[i] for i in range(len(r0))]
        r4_m = [(T2_m[i] - 5*r6_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r4_m[i] = (r4_m[i] - 21*r8_m[i]) % m
        
        r2 = [T1[i] - r4[i] - r6[i] - r8[i] for i in range(len(r0))]
        r2_m = [(T1_m[i] - r4_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r2_m[i] = (r2_m[i] - r6_m[i]) % m
            r2_m[i] = (r2_m[i] - r8_m[i]) % m
    
        T7 = [(r[1][i] - r[-1][i])//2 for i in range(len(r0))]
        T7_m = [(r_m[1][i] - r_m[-1][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T7_m[i] = T7_m[i] // 2
    
        T8 = [((r[2][i] - r[-2][i])//4
              - T7[i]) // 3 for i in range(len(r0))]
        T8_m = [(r_m[2][i] - r_m[-2][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T8_m[i] = T8_m[i] // 4
            T8_m[i] = (T8_m[i] - T7_m[i]) % m
            T8_m[i] = (T8_m[i] * inv3) % m
    
        T9 = [((r[3][i] - r[-3][i])//6
              - T7[i]) // 8 for i in range(len(r0))]
        T9_m = [(r_m[3][i] - r_m[-3][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T9_m[i] = T9_m[i] // 2
            T9_m[i] = (T9_m[i] * inv3) % m
            T9_m[i] = (T9_m[i] - T7_m[i]) % m
            T9_m[i] = T9_m[i] // 8
            
        T10 = [((r[4][i] - r[-4][i])//8
              - T7[i]) // 15 for i in range(len(r0))]
        T10_m = [(r_m[4][i] - r_m[-4][i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T10_m[i] = T10_m[i] // 8
            T10_m[i] = (T10_m[i] - T7_m[i]) % m
            T10_m[i] = (T10_m[i] * inv15) % m
    
        T11 = [((r[5][i] - r0[i] - 25*r2[i] - 625*r4[i] 
                       - 15625*r6[i] - 390625*r8[i] - 9765625*r10[i])//5
             - T7[i]) // 24 for i in range(len(r0))]
        T11_m = [(r_m[5][i] - r0_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T11_m[i] = (T11_m[i] - 25*r2_m[i]) % m
            T11_m[i] = (T11_m[i] - 625*r4_m[i]) % m
            T11_m[i] = (T11_m[i] - 15625*r6_m[i]) % m
            T11_m[i] = (T11_m[i] - 390625*r8_m[i]) % m
            T11_m[i] = (T11_m[i] - 9765625*r10_m[i]) % m
            T11_m[i] = (T11_m[i] * inv5) % m
            T11_m[i] = (T11_m[i] - T7_m[i]) % m
            T11_m[i] = T11_m[i] // 8
            T11_m[i] = (T11_m[i] * inv3) % m
    
        T12 = [(T9[i] - T8[i]) // 5 for i in range(len(r0))]
        T12_m = [((T9_m[i] - T8_m[i]) * inv5) % m for i in range(len(r0_m))]
    
        T13 = [(T10[i] - T9[i]) // 7 for i in range(len(r0))]
        T13_m = [((T10_m[i] - T9_m[i]) * inv7) % m for i in range(len(r0_m))]
        
        T14 = [(T11[i] - T10[i]) // 9 for i in range(len(r0))]
        T14_m = [((T11_m[i] - T10_m[i]) * inv9) % m for i in range(len(r0_m))]
        
        T15 = [(T13[i] - T12[i]) // 12 for i in range(len(r0))]
        T15_m = [(T13_m[i] - T12_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T15_m[i] = T15_m[i] // 4
            T15_m[i] = (T15_m[i] * inv3) % m
            
        T16 = [(T14[i] - T13[i]) // 16 for i in range(len(r0))]
        T16_m = [(T14_m[i] - T13_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            T16_m[i] = T16_m[i] // 16
            
        r9 = [(T16[i] - T15[i]) // 21 for i in range(len(r0))]
        r9_m = [(T16_m[i] - T15_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r9_m[i] = (r9_m[i] * inv21) % m
            
        r7 = [T15[i] - 30*r9[i] for i in range(len(r0))]
        r7_m = [(T15_m[i] - 30*r9_m[i]) % m for i in range(len(r0_m))]            
    
        r5 = [T12[i] - 14*r7[i] - 147*r9[i] for i in range(len(r0))]
        r5_m = [(T12_m[i] - 14*r7_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r5_m[i] = (r5_m[i] - 147*r9_m[i]) % m
    
        r3 = [T8[i] - 5*r5[i] - 21*r7[i] - 85*r9[i] for i in range(len(r0))]
        r3_m = [(T8_m[i] - 5*r5_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r3_m[i] = (r3_m[i] - 21*r7_m[i]) % m
            r3_m[i] = (r3_m[i] - 85*r9_m[i]) % m
    
        r1 = [T7[i] - r3[i] - r5[i] - r7[i] - r9[i] for i in range(len(r0))]
        r1_m = [(T7_m[i] - r3_m[i]) % m for i in range(len(r0_m))]
        for i in range(len(r0_m)):
            r1_m[i] = (r1_m[i] - r5_m[i]) % m
            r1_m[i] = (r1_m[i] - r7_m[i]) % m
            r1_m[i] = (r1_m[i] - r9_m[i]) % m
    
        r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
        r_coefs_m = (r0_m, r1_m, r2_m, r3_m, r4_m, r5_m, r6_m, r7_m, r8_m,
                     r9_m, r10_m)
        
        T_coefs = (T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13,
                   T14, T15, T16)
        T_coefs_m = (T1_m, T2_m, T3_m, T4_m, T5_m, T6_m, T7_m, T8_m, T9_m,
                     T10_m, T11_m, T12_m, T13_m, T14_m, T15_m, T16_m)
        
        for i in eval_list:
            label_string = "r({})".format(i)
            D[label_string].append(exp - max_power_of_two_list(r[i], r_m[i], exp))
        for i in range(2*6 - 1):
            label_string = "r{}".format(i)
            D[label_string].append(exp - max_power_of_two_list(r_coefs[i], r_coefs_m[i], exp))
        for i in range(1, 17):
             label_string = "T{}".format(i)
             D[label_string].append(exp - max_power_of_two_list(T_coefs[i-1], T_coefs_m[i-1], exp))
    for i in eval_list:
        label_string = "r({})".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(2*6 - 1):
        label_string = "r{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    for i in range(1,17):
        label_string = "T{}".format(i)
        print(label_string + " : {} ".format(stats.mode(D[label_string])[0][0]))
    return            
            
check_toom4_cost(num_trials=50)           
            