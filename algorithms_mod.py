# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:43:52 2019

@author: Marcus
Now that I know what I'm doing, this does everything mod 2048.
"""

# need to take a ceiling, that's all
import numpy as np

def multiply_2048(f, g, algorithm, minsize=1, using_toom_6=False, outer=True):
    """ This multiplies f and g. algorithm_list specifies how.  Our valid
        algorithms are 2, 3, 4, 5, and 6, corresponding to Toom-i (viewing
        Karatsuba as Toom-2). Update: algorithm must be a list. 
        If algorithm_list is a list of valid ints, then we recursively 
        multiply until the end of the list, at which point we use 
        schoolbook.  Also, 1 means schoolbook.
        This does it mod 2048. All integers in the lists should be
        between 0 and 2047. Each intermediate step is mod 66536,
        unless we use Toom-6. Then, we need to use bigger ints.
        The parameter outer is for knowing if this is the first algorithm in
        a list of algorithms. If it is, then the final mod is mod 2048."""
       
    if len(f) != len(g):
        raise ValueError("Can only multiply polys of the same length")
        
    if algorithm == 6 or (type(algorithm) == list and 6 in algorithm):
        using_toom_6 = True
        
    if using_toom_6:  # in order to divide by 128, we need 18 bit ints
        m = 2**(32)   # to simulate C++, we are using 32 bit ints
    else:
        m = 2**(16)   # otherwise, we can use 16-bit ints.
    
    # do schoolbook multiplication if the list is empty, or
    # if the polynomials are shorter than the Toom specified
    if ((algorithm in range(2,7) and 
         (len(f) < algorithm or len(f) < minsize)) 
        or 
        (type(algorithm) == list 
         and len(algorithm) == 0)
        or algorithm == 1):
            if not outer:
                return schoolbook_mod(f, g, m)
            else:
                return schoolbook_mod(f, g, 2048)
    
    if type(algorithm) == list:
        n = algorithm[0]
        next_alg = algorithm[1:]
    else:
        n = algorithm
        next_alg = algorithm
        
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_eval_list(n)
    
    # plug the numbers in
    f_eval = evaluate_blocks_list_mod(fblocks, eval_list, m)
    g_eval = evaluate_blocks_list_mod(gblocks, eval_list, m)
    
    # perform the recursive multiplication
    r = {eval_list[i]:multiply_2048(f_eval[i], 
                                    g_eval[i], 
                                    next_alg, 
                                    minsize=minsize,
                                    using_toom_6=using_toom_6,
                                    outer=False)
            for i in range(len(f_eval))}
    #print(r)
    #for i in r:
        #print(i, [x % 65536 for x in r[i]])
    # Solve for the coefficients
    r_coefs = solve_for_coefficients_mod(n, r, m)
    
    #for i in range(len(r_coefs)):
        #print(i, [x % 65536 for x in r_coefs[i]])
        #print(i, [x % 2048 for x in r_coefs[i]])
    
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
    
    #print(f, " times ", g, " = ", [x % 65536 for x in prod[:2*len(f)-1]])

    if outer:
        return [x % 2048 for x in prod[:2*len(f)-1]]
    else:
        return [x % m for x in prod[:2*len(f)-1]]
    


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


def solve_for_coefficients_mod(n, r, m):
    """ This function handles the long formulas needed to explicitly find
        the Toom formulas for Toom-2 up to Toom-6."""
    
    # manually compute inverses, since these are the only ones we
    # will need
    if m == (2**16):
        inv3 = 43691
        inv15 = 61167
        inv45 = 20389
        inv315 = 12275
        inv945 = 25937
        inv14175 = 58527
    elif m == (2**32):
        inv3 = 2863311531
        inv15 = 4008636143
        inv45 = 2767867813
        inv315 = 4076810227
        inv945 = 4222248273
        inv14175 = 3717457055
    else:
        raise ValueError("what is m")
    
    if n == 2:
                
        r0 = r[0]
        
        r2 = r['infinity']
        
        r1 = [r[1][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r1[i] = (r1[i] - r0[i]) % m
            r1[i] = (r1[i] - r2[i]) % m
            
        return (r0, r1, r2)
    
    if n == 3:
               
        r0 = r[0]
        
        r4 = r['infinity']
        
        r2 = [(r[1][i] + r[-1][i])//2 for i in range(len(r0))]
        for i in range(len(r0)):
            r2[i] = (r2[i] - r0[i]) % m
            r2[i] = (r2[i] - r4[i]) % m
        r3 = [r[2][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r3[i] = (r3[i] - r0[i]) % m
            r3[i] = (r3[i] - 4*r2[i]) % m
            r3[i] = (r3[i] - 16*r4[i]) % m
            r3[i] = (r3[i] + r[-1][i]) % m
            r3[i] = (r3[i] - r[1][i]) % m
            r3[i] = r3[i]//2
            r3[i] = (r3[i] * inv3) % m
        
        r1 = [r[1][i] - r0[i] for i in range(len(r0))]
        for i in range(len(r0)):
            r1[i] = (r1[i] - r2[i]) % m
            r1[i] = (r1[i] - r3[i]) % m
            r1[i] = (r1[i] - r4[i]) % m
        
        return (r0, r1, r2, r3, r4)
    
    if n == 4:
        
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
      
        print([x % 2048 for x in r1])
        return (r0, r1, r2, r3, r4, r5, r6)
    
    if n == 5:
        
        r0 = r[0]
    
        r8 = r['infinity']
    
        r6 = [r[3][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r6[i] = (r6[i] + r[-3][i]) % m
            r6[i] = (r6[i] - 6*r[2][i]) % m
            r6[i] = (r6[i] - 6*r[-2][i]) % m
            r6[i] = (r6[i] + 15*r[1][i]) % m
            r6[i] = (r6[i] + 15*r[-1][i]) % m
            r6[i] = (r6[i] - 20*r0[i]) % m
            r6[i] = (r6[i] - 10080*r8[i]) % m
            r6[i] = r6[i] // 16
            r6[i] = (r6[i] * inv45) % m
      
        r4 = [r[2][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r4[i] = (r4[i] + r[-2][i]) % m
            r4[i] = (r4[i] - 4*r[1][i]) % m
            r4[i] = (r4[i] - 4*r[-1][i]) % m
            r4[i] = (r4[i] + 6*r0[i]) % m
            r4[i] = (r4[i] - 120*r6[i]) % m
            r4[i] = (r4[i] - 504*r8[i]) % m
            r4[i] = r4[i] // 8
            r4[i] = (r4[i] * inv3) % m
        
        r2 = [r[1][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r2[i] = (r2[i] + r[-1][i]) % m
            r2[i] = (r2[i] - 2*r0[i]) % m
            r2[i] = (r2[i] - 2*r4[i]) % m
            r2[i] = (r2[i] - 2*r6[i]) % m
            r2[i] = (r2[i] - 2*r8[i]) % m
            r2[i] = r2[i] // 2
        
        r7 = [r[4][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r7[i] = (r7[i] - 14*r[1][i]) % m
            r7[i] = (r7[i] + 14*r[2][i]) % m
            r7[i] = (r7[i] - 6*r[3][i]) % m
            r7[i] = (r7[i] + 5*r0[i]) % m
            r7[i] = (r7[i] - 4*r2[i]) % m
            r7[i] = (r7[i] + 20*r4[i]) % m
            r7[i] = (r7[i] - 604*r6[i]) % m
            r7[i] = (r7[i] - 29740*r8[i]) % m
            r7[i] = r7[i] // 16
            r7[i] = (r7[i] * inv315) % m
        
        r5 = [r[3][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r5[i] = (r5[i] + 5*r[1][i]) % m
            r5[i] = (r5[i] - 4*r[2][i]) % m
            r5[i] = (r5[i] - 2*r0[i]) % m
            r5[i] = (r5[i] + 2*r2[i]) % m
            r5[i] = (r5[i] - 22*r4[i]) % m
            r5[i] = (r5[i] - 478*r6[i]) % m
            r5[i] = (r5[i] - 1680*r7[i]) % m
            r5[i] = (r5[i] - 5542*r8[i]) % m
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
            r3[i] = (r3[i] - 126*r7[i]) % m
            r3[i] = (r3[i] - 254*r8[i]) % m
            r3[i] = r3[i] // 2
            r3[i] = (r3[i] * inv3) % m
    
        r1 = [r[1][i] - r0[i] for i in range(len(r0))]
        for i in range(len(r0)):
            r1[i] = (r1[i] - r2[i]) % m
            r1[i] = (r1[i] - r3[i]) % m
            r1[i] = (r1[i] - r4[i]) % m
            r1[i] = (r1[i] - r5[i]) % m
            r1[i] = (r1[i] - r6[i]) % m
            r1[i] = (r1[i] - r7[i]) % m
            r1[i] = (r1[i] - r8[i]) % m
       
        return (r0, r1, r2, r3, r4, r5, r6, r7, r8)
            
    if n == 6:
        
        r0 = r[0]
    
        r10 = r['infinity']
    
        r8 = [r[4][i] + r[-4][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r8[i] = (r8[i] - 8*r[3][i]) % m
            r8[i] = (r8[i] - 8*r[-3][i]) % m
            r8[i] = (r8[i] + 28*r[2][i]) % m
            r8[i] = (r8[i] + 28*r[-2][i]) % m
            r8[i] = (r8[i] - 56*r[1][i]) % m
            r8[i] = (r8[i] - 56*r[-1][i]) % m
            r8[i] = (r8[i] + 70*r0[i]) % m
            r8[i] = (r8[i] - 1209600*r10[i]) % m
            r8[i] = r8[i] // 128
            r8[i] = (r8[i] * inv315) % m
    
        r6 = [r[3][i] + r[-3][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r6[i] = (r6[i] - 6*r[2][i]) % m
            r6[i] = (r6[i] - 6*r[-2][i]) % m
            r6[i] = (r6[i] + 15*r[1][i]) % m
            r6[i] = (r6[i] + 15*r[-1][i]) % m
            r6[i] = (r6[i] - 20*r0[i]) % m
            r6[i] = (r6[i] - 10080*r8[i]) % m
            r6[i] = (r6[i] - 105840*r10[i]) % m
            r6[i] = r6[i] // 16
            r6[i] = (r6[i] * inv45) % m
        
        r4 = [r[2][i] + r[-2][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r4[i] = (r4[i] - 4*r[1][i]) % m
            r4[i] = (r4[i] - 4*r[-1][i]) % m
            r4[i] = (r4[i] + 6*r0[i]) % m
            r4[i] = (r4[i] - 120*r6[i]) % m
            r4[i] = (r4[i] - 504*r8[i]) % m
            r4[i] = (r4[i] - 2040*r10[i]) % m
            r4[i] = r4[i] // 8
            r4[i] = (r4[i] * inv3) % m
        
        r2 = [(r[1][i] + r[-1][i]) // 2 for i in range(len(r0))]
        for i in range(len(r0)):
            r2[i] = (r2[i] - r0[i]) % m
            r2[i] = (r2[i] - r4[i]) % m
            r2[i] = (r2[i] - r6[i]) % m
            r2[i] = (r2[i] - r8[i]) % m
            r2[i] = (r2[i] - r10[i]) % m
            
        r9 = [(5*r[5][i]) % m for i in range(len(r0))]    
        for i in range(len(r0)):
            r9[i] = (r9[i] - 40*r[4][i]) % m
            r9[i] = (r9[i] + 135*r[3][i]) % m
            r9[i] = (r9[i] - 240*r[2][i]) % m
            r9[i] = (r9[i] + 210*r[1][i]) % m
            r9[i] = (r9[i] - 70*r0[i]) % m
            r9[i] = (r9[i] + 50*r2[i]) % m
            r9[i] = (r9[i] - 190*r4[i]) % m
            r9[i] = (r9[i] + 2450*r6[i]) % m
            r9[i] = (r9[i] - 156190*r8[i]) % m
            r9[i] = (r9[i] - 14611150*r10[i]) % m
            r9[i] = r9[i] // 128
            r9[i] = (r9[i] * inv14175) % m
        
        r7 = [r[5][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r7[i] = (r7[i] - 2*r[4][i]) % m
            r7[i] = (r7[i] - 9*r[3][i]) % m
            r7[i] = (r7[i] + 36*r[2][i]) % m
            r7[i] = (r7[i] - 42*r[1][i]) % m
            r7[i] = (r7[i] + 16*r0[i]) % m
            r7[i] = (r7[i] - 14*r2[i]) % m
            r7[i] = (r7[i] + 82*r4[i]) % m
            r7[i] = (r7[i] - 3134*r6[i]) % m
            r7[i] = (r7[i] - 209678*r8[i]) % m
            r7[i] = (r7[i] - 1270080*r9[i]) % m
            r7[i] = (r7[i] - 7173854*r10[i]) % m
            r7[i] = r7[i] // 32
            r7[i] = (r7[i] * inv945) % m
        
        r5 = [r[3][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r5[i] = (r5[i] - 4*r[2][i]) % m
            r5[i] = (r5[i] + 5*r[1][i]) % m
            r5[i] = (r5[i] - 2*r0[i]) % m
            r5[i] = (r5[i] + 2*r2[i]) % m
            r5[i] = (r5[i] - 22*r4[i]) % m
            r5[i] = (r5[i] - 478*r6[i]) % m
            r5[i] = (r5[i] - 1680*r7[i]) % m
            r5[i] = (r5[i] - 5542*r8[i]) % m
            r5[i] = (r5[i] - 17640*r9[i]) % m
            r5[i] = (r5[i] - 54958*r10[i]) % m
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
            r3[i] = (r3[i] - 126*r7[i]) % m
            r3[i] = (r3[i] - 254*r8[i]) % m
            r3[i] = (r3[i] - 510*r9[i]) % m
            r3[i] = (r3[i] - 1022*r10[i]) % m
            r3[i] = r3[i] // 2
            r3[i] = (r3[i] * inv3) % m
        
        r1 = [r[1][i] for i in range(len(r0))]
        for i in range(len(r0)):
            r1[i] = (r1[i] - r0[i]) % m
            r1[i] = (r1[i] - r2[i]) % m
            r1[i] = (r1[i] - r3[i]) % m
            r1[i] = (r1[i] - r4[i]) % m
            r1[i] = (r1[i] - r5[i]) % m
            r1[i] = (r1[i] - r6[i]) % m
            r1[i] = (r1[i] - r7[i]) % m
            r1[i] = (r1[i] - r8[i]) % m
            r1[i] = (r1[i] - r9[i]) % m
            r1[i] = (r1[i] - r10[i]) % m
                
        return (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
     
    else:
        raise ValueError("We haven't implemented Toom-" + str(n))
        
def number_of_twos_in_n(n):
    """ Returns the highest power of 2 dividing n"""
    exp = 0
    power = 2
    while n % power == 0:
        exp += 1
        power *= 2
    return exp

def number_of_twos_in_factorial(n):
    """ Returns the highest power of 2 dividing n! """
    num = 0
    for i in range(1, n + 1):
        num += number_of_twos_in_n(i)
    return num

    
def algorithm_cost(algorithm):
    """ algorithm is a list specifying the Toom algorithm. This counts
        the total number of powers of two that we lose by using this
        algorithm list"""
    cost = 0
    for n in algorithm:
        # in Toom-2 we don't have to divide
        if n != 2:
            cost += number_of_twos_in_factorial(2*n - 3)
    return cost

def max_power_of_two(a, b, max_n=11):
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

print(algorithm_cost([12]))

f = range(1, 745)
g = [745 - x for x in range(1,745)]
X = multiply_2048(f, g, [4, 3, 2])
#A = multiply_2048(f, g, [6])
#B = schoolbook_mod(f, g, 2048)
#print(max_power_of_two_list(A,B, max_n=32))