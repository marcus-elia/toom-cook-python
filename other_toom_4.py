# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:20:37 2019

@author: Marcus
This is for seeing the cost of the Toom-4 interpolation
formulas described by Schwabe et al
"""

import numpy as np

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

def make_eval_list(n):
    """ In Toom-n, this makes the list of number to plug in"""
    if n == 2:
        return (0, 1, 'infinity')
    if n == 3:
        return (0, 1, -1, 2, 'infinity')
    if n == 4:
        return (0, 1, -1, 2, -2, 3, 'infinity')
    if n == 5:
        return (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    if n == 6:
        return (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    else:
        raise ValueError("We haven't implemented Toom-" + str(n))
        
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

def plain_toom_4_mod(f, g, m):
    n = 4
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list_mod(fblocks, eval_list, m)
    g_eval = evaluate_blocks_list_mod(gblocks, eval_list, m)

    r = {eval_list[i] : schoolbook_mod(f_eval[i], g_eval[i], m)
                for i in range(len(f_eval))}
    
    
    inv3 = 43691
    inv15 = 61167
    
    r0 = r[0]
        
    r6 = r['infinity']
     
    r4 = [6*r0[i] % m for i in range(len(r0))]
    for i in range(len(r0)):
        r4[i] = (r4[i] - 120*r6[i]) % m
        r4[i] = (r4[i] + r[2][i] + r[-2][i]) % m
        r4[i] = (r4[i] - 4*(r[1][i] + r[-1][i])) % m
        r4[i] = r4[i] // 8
        r4[i] = (r4[i] * inv3) % m
                        
    r2 = [r[1][i] + r[-1][i] for i in range(len(r0))]
    for i in range(len(r0)):
        r2[i] = (r2[i] - 2*r0[i]) % m
        r2[i] = (r2[i] - 2*r4[i]) % m
        r2[i] = (r2[i] - 2*r6[i]) % m
        r2[i] = r2[i] // 2
            
    r5 = [r[3][i] for i in range(len(r0))]
    for i in range(len(r0)):
        r5[i] = (r5[i] - 4*r[2][i]) % m
        r5[i] = (r5[i] + 5*r[1][i]) % m
        r5[i] = (r5[i] - 2*r0[i]) % m
        r5[i] = (r5[i] + 2*r2[i]) % m
        r5[i] = (r5[i] - 22*r4[i]) % m
        r5[i] = (r5[i] - 478*r6[i]) % m
        r5[i] = r5[i] // 8
        r5[i] = (r5[i] * inv15) % m
        
    r3 = [r[2][i] for i in range(len(r0))]
    for i in range(len(r0)):
        r3[i] = (r3[i] - 2*r[1][i]) % m
        r3[i] = (r3[i] + r0[i]) % m
        r3[i] = (r3[i] - 2*r2[i]) % m
        r3[i] = (r3[i] - 14*r4[i]) % m
        r3[i] = (r3[i] - 30*r5[i]) % m
        r3[i] = (r3[i] - 62*r6[i]) % m
        r3[i] = r3[i] // 2
        r3[i] = (r3[i] * inv3) % m
     
    r1 = [r[1][i] - r0[i] for i in range(len(r0))]
    for i in range(len(r0)):
        r1[i] = (r1[i] - r2[i]) % m
        r1[i] = (r1[i] - r3[i]) % m
        r1[i] = (r1[i] - r4[i]) % m
        r1[i] = (r1[i] - r5[i]) % m
        r1[i] = (r1[i] - r6[i]) % m
        
    r_coefs = (r0, r1, r2, r3, r4, r5, r6)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [(r_coefs[j-1][k+i] + r_coefs[j][i])
                       for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [(r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i])
                   for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return [x % m for x in prod[:2*len(f)-1]]

def toom_4(f, g):
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
    #V_1 = [(r[2][i] - r[-2][i]) for i in range(len(r0))]
    #for i in range(len(r0)):
    #    V_1[i] = V_1[i] // 4 
    #    V_1[i] = (V_1[i] - T2[i])
    #    V_1[i] = (V_1[i] // 3)
    #V_2 = [(r[3][i] - r0[i] 
    #               - 9*r2[i] 
     #              - 81*r4[i] 
     #              - 729*r6[i]) for i in range(len(r0))]
    #V_2 = [(V_2[i] - T2[i]) for i in range(len(r0))]
    #for i in range(len(r0)):
    #    V_2[i] = V_2[i] // 8
    #    V_2[i] = (V_2[i] - V_1[i])
    #print(V_2)
    T3 = [((r[2][i] - r[-2][i]) // 4 - T2[i]) // 3 for i in L]
    
    T4 = [((r[3][i] - r0[i] - 9*r2[i] - 81*r4[i] - 729*r6[i]) // 3
           - T2[i]) // 8
         for i in L]

    r5 = [(T4[i] - T3[i]) // 5 for i in L]
    
    r3 = [T3[i] - 5*r5[i] for i in L]
    
    r1 = [T2[i] - r3[i] - r5[i] for i in L]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6)
    for h in r_coefs:
       print(h)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]

def other_toom_4(f, g):
    """ m is what we are working mod (like 2^16 for 16-bit ints)"""
    n = 4
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_eval_list(n)
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    C = {eval_list[i]:schoolbook(f_eval[i], 
                                    g_eval[i])
            for i in range(len(f_eval))}
    
    C_0 = C[0]
    
    C_6 = C['infinity']
    
    V_0 = [(C[1][i] + C[-1][i]) // 2 for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_0[i] = (V_0[i] - C_0[i])
        V_0[i] = (V_0[i] - C_6[i])
    
    V_1 = [(C[2][i] + C[-2][i] - 2*C_0[i]- 128*C_6[i]) for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_1[i] = V_1[i] // 8
    
    C_4 = [(V_1[i] - V_0[i]) for i in range(len(C_0))]
    for i in range(len(C_0)):
        C_4[i] = (C_4[i] // 3)
    
    C_2 = [(V_0[i] - C_4[i]) for i in range(len(C_0))]
    
    V_0 = [(C[1][i] - C[-1][i]) for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_0[i] = V_0[i] // 2
        
    V_1 = [(C[2][i] - C[-2][i]) for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_1[i] = V_1[i] // 4 
        V_1[i] = (V_1[i] - V_0[i])
        V_1[i] = (V_1[i] // 3)

    V_2 = [(C[3][i] - C_0[i] 
                   - 9*C_2[i] 
                   - 81*C_4[i] 
                   - 729*C_6[i]) for i in range(len(C_0))]
   
    V_2 = [(V_2[i] - 3*V_0[i]) for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_2[i] = V_2[i] // 8
        V_2[i] = (V_2[i] - 3*V_1[i])
        V_2[i] = V_2[i] // 3

    C_5 = [(V_2[i] // 5) for i in range(len(C_0))]
    
    C_3 = [(V_1[i] - V_2[i]) for i in range(len(C_0))]
    
    C_1 = [(V_0[i] - C_3[i] - C_5[i]) for i in range(len(C_0))]
    
    C_coefs = (C_0, C_1, C_2, C_3, C_4, C_5, C_6)
    for h in C_coefs:
        print(h)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = C_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [(C_coefs[j-1][k+i] + C_coefs[j][i])
                       for i in range(k-1)]
        prod = prod + [C_coefs[j][k-1]]

    prod = prod + [(C_coefs[2*n-3][k+i] + C_coefs[2*n-2][i])
                   for i in range(k-1)]

    prod = prod + C_coefs[2*n-2][k-1:]

    return [x for x in prod[:2*len(f)-1]]

def other_toom_4_mod(f, g, m):
    """ m is what we are working mod (like 2^16 for 16-bit ints)"""
    n = 4
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_eval_list(n)
    
    # plug the numbers in
    f_eval = evaluate_blocks_list_mod(fblocks, eval_list, m)
    g_eval = evaluate_blocks_list_mod(gblocks, eval_list, m)
    
    # perform the recursive multiplication
    C = {eval_list[i]:schoolbook_mod(f_eval[i], 
                                    g_eval[i], 
                                    m)
            for i in range(len(f_eval))}
    
    inv3 = 43691
    inv5 = 52429
    
    C_0 = C[0]
    
    C_6 = C['infinity']
    
    V_0 = [(C[1][i] + C[-1][i]) // 2 for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_0[i] = (V_0[i] - C_0[i]) % m
        V_0[i] = (V_0[i] - C_6[i]) % m
    
    V_1 = [(C[2][i] + C[-2][i] - 2*C_0[i]- 128*C_6[i]) % m for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_1[i] = V_1[i] // 8
    
    C_4 = [(V_1[i] - V_0[i]) % m for i in range(len(C_0))]
    for i in range(len(C_0)):
        C_4[i] = (C_4[i] * inv3) % m
    
    C_2 = [(V_0[i] - C_4[i]) % m for i in range(len(C_0))]
    
    V_0 = [(C[1][i] - C[-1][i]) % m for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_0[i] = V_0[i] // 2
    
    V_1 = [(C[2][i] - C[-2][i]) % m  for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_1[i] = V_1[i] // 4 
        V_1[i] = (V_1[i] - V_0[i]) % m
        V_1[i] = (V_1[i] * inv3) % m
    
    V_2 = [(C[3][i] - C_0[i] 
                   - 9*C_2[i] 
                   - 81*C_4[i] 
                   - 729*C_6[i]) % m for i in range(len(C_0))]
    
    V_2 = [(V_2[i] - 3*V_0[i]) % m for i in range(len(C_0))]
    for i in range(len(C_0)):
        V_2[i] = V_2[i] // 8
        V_2[i] = (V_2[i] - 3*V_1[i]) % m
        V_2[i] = (V_2[i] * inv3) % m
        
    C_5 = [(V_2[i] * inv5) % m for i in range(len(C_0))]
    
    C_3 = [(V_1[i] - V_2[i]) % m for i in range(len(C_0))]
    
    C_1 = [(V_0[i] - C_3[i] - C_5[i]) % m for i in range(len(C_0))]
    
    C_coefs = (C_0, C_1, C_2, C_3, C_4, C_5, C_6)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = C_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [(C_coefs[j-1][k+i] + C_coefs[j][i])
                       for i in range(k-1)]
        prod = prod + [C_coefs[j][k-1]]

    prod = prod + [(C_coefs[2*n-3][k+i] + C_coefs[2*n-2][i])
                   for i in range(k-1)]

    prod = prod + C_coefs[2*n-2][k-1:]

    return [x % m for x in prod[:2*len(f)-1]]


#f = [int(x) for x in np.random.randint(0, 2**(12), size=4)]
#g = [int(x) for x in np.random.randint(0, 2**(12), size=4)]
f = [1,2,3,4]
g = [4,3,2,1]
good = schoolbook_mod(f, g, 2**(16))
fake = other_toom_4_mod(f, g, 2**(16))
print(good)
print([f % (2**13) for f in fake])
print([(good[i] - fake[i]) % (2**(13)) for i in range(len(good))])