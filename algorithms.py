# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 14:00:04 2019

@author: Marcus

This is my most pristine version of combining different
Toom-Cook algorithms to multiply polynomials.
"""

# need to take a ceiling, that's all
import numpy as np

def multiply(f, g, algorithm):
    """ This multiplies f and g. algorithm_list specifies how.  Our valid
        algorithms are 2, 3, 4, 5, 6, 8 corresponding to Toom-i (viewing
        Karatsuba as Toom-2). Update: algorithm must be a list.
        If algorithm_list is a list of valid ints, then we recursively 
        multiply until the end of the list, at which point we use 
        schoolbook.  Also, 1 means schoolbook."""
       
    if len(f) != len(g):
        raise ValueError("Can only multiply polys of the same length")
        
    # do schoolbook multiplication if the list is empty
    if len(algorithm) == 0:
        return schoolbook(f, g)
    
    
    n = algorithm[0]
    next_alg = algorithm[1:]
            
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_eval_list(n)
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    r = {eval_list[i]:multiply(f_eval[i], g_eval[i], next_alg)
            for i in range(len(f_eval))}
    
    # Solve for the coefficients    
    r_coefs = solve_for_coefficients(n, r)
    
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]
    


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
    """ In Toom-n, this makes the list of number to plug in
        (0, 1, -1, 2, -2, ..., n-2, -(n-2), n-1, 'infinity')"""
    eval_points = [0]
    for i in range(1, n - 1):
        eval_points.append(i)
        eval_points.append(-i)
    eval_points.append(n-1)
    eval_points.append('infinity')
    return eval_points
  

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


def solve_for_coefficients(n, r):
    """ This function handles the long formulas needed to explicitly find
        the Toom formulas for Toom-2 up to Toom-6."""
    if n == 2:
        r0 = r[0]
        r2 = r['infinity']
        r1 = [r[1][i] - r0[i] - r2[i] for i in range(len(r0))]
        return (r0, r1, r2)
    
    if n == 3:
        
        r0 = r[0]
        
        r4 = r['infinity']
        
        r2 = [(r[1][i] + r[-1][i])//2 
              - r0[i] 
              - r4[i] 
              for i in range(len(r0))]
        print([x for x in r[1]])
        print([x for x in r[-1]])
        print([((r[1][i] + r[-1][i])//2) % (2**(16)) for i in range(len(r[1]))])
        r3 = [(r[2][i] 
              - r0[i] 
              - 4*r2[i] 
              - 16*r4[i] 
              + r[-1][i] 
              - r[1][i])//6 
             for i in range(len(r0))]
        
        r1 = [r[1][i] - r0[i] - r2[i] - r3[i] - r4[i] for i in range(len(r0))]
        
        return (r0, r1, r2, r3, r4)
    
    if n == 4:
        
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
        
        return (r0, r1, r2, r3, r4, r5, r6)
    
    if n == 5:
        
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
               - 126*r7[i] 
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
    
        return (r0, r1, r2, r3, r4, r5, r6, r7, r8)
     
    if n == 6:
        
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
    
        #r9 = [(5*r[5][i]
        #       - 40*r[4][i] 
        #       + 135*r[3][i] 
        #       - 240*r[2][i] 
        #       + 210*r[1][i]
        #       - 70*r0[i] 
        #       + 50*r2[i] 
        #       - 190*r4[i] 
        #       + 2450*r6[i] 
        #       - 156190*r8[i] 
        #       - 14611150*r10[i])//1814400 
        #     for i in range(len(r0))]
        
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
            ) // 362880 for i in range(len(r0))]
    
        #r7 = [(r[5][i] 
        #       - 2*r[4][i] 
        #       - 9*r[3][i] 
        #       + 36*r[2][i] 
        #       - 42*r[1][i] 
        #       + 16*r0[i] 
        #       - 14*r2[i] 
        #       + 82*r4[i] 
        #       - 3134*r6[i] 
        #       - 209678*r8[i] 
        #       - 1270080*r9[i] 
        #       - 7173854*r10[i])//30240 
        #      for i in range(len(r0))]
    
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
              ) // 5040 for i in range(len(r0))]
    
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
        
        return (r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
     
    if n == 8:
        
        r0 = r[0]

        r14 = r['infinity']

        r12 = [(-792*(r[1][i] + r[-1][i])
                + 495*(r[2][i] + r[-2][i])
                - 220*(r[3][i] + r[-3][i])
                + 66*(r[4][i] + r[-4][i])
                - 12*(r[5][i] + r[-5][i])
                + (r[6][i] + r[-6][i])
                + 924*r0[i]
                - 43589145600*r14[i]
               ) // 479001600 for i in range(len(r0))]

        r10 = [(210*(r[1][i] + r[-1][i])
                - 120*(r[2][i] + r[-2][i])
                + 45*(r[3][i] + r[-3][i])
                - 10*(r[4][i] + r[-4][i])
                + (r[5][i] + r[-5][i])
                - 252*r0[i]
                - 199584000*r12[i]
                - 7264857600*r14[i]
               ) // 3628800 for i in range(len(r0))]

        r8 = [(-56*(r[1][i] + r[-1][i])
                + 28*(r[2][i] + r[-2][i])
                - 8*(r[3][i] + r[-3][i])
                + (r[4][i] + r[-4][i])
                + 70*r0[i]
                - 1209600*r10[i]
                - 25280640*r12[i]
                - 461260800*r14[i]
               ) // 40320 for i in range(len(r0))]

        r6 = [(15*(r[1][i] + r[-1][i])
                - 6*(r[2][i] + r[-2][i])
                + (r[3][i] + r[-3][i])
                - 20*r0[i]
                - 10080*r8[i]
                - 105840*r10[i]
                - 1013760*r12[i]
                - 9369360*r14[i]
               ) // 720 for i in range(len(r0))]

        r4 = [(-4*(r[1][i] + r[-1][i])
                + (r[2][i] + r[-2][i])
                + 6*r0[i]
                - 120*r6[i]
                - 504*r8[i]
                - 2040*r10[i]
                - 8184*r12[i]
                - 32760*r14[i]
               ) // 24 for i in range(len(r0))]

        r2 = [((r[1][i] + r[-1][i])
                - 2*r0[i]
                - 2*r4[i]
                - 2*r6[i]
                - 2*r8[i]
                - 2*r10[i]
                - 2*r12[i]
                - 2*r14[i]
               ) // 2 for i in range(len(r0))]

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
                 for i in range(len(r0))]
        
        return (r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14)
        
    if n == 12:
        
        r0 = r[0]

        r22 = r['infinity']

        r20 = [(-167960*(r[1][i] + r[-1][i])
               + 125970*(r[2][i] + r[-2][i])
               - 77520*(r[3][i] + r[-3][i])
               + 38760*(r[4][i] + r[-4][i])
               - 15504*(r[5][i] + r[-5][i])
               + 4845*(r[6][i] + r[-6][i])
               - 1140*(r[7][i] + r[-7][i])
               + 190*(r[8][i] + r[-8][i])
               - 20*(r[9][i] + r[-9][i])
               + (r[10][i] + r[-10][i])
               + 184756*r0[i]
               - 936667273148006400000*r22[i]
              ) // 2432902008176640000 for i in range(len(r0))]

        r18 = [(43758*(r[1][i] + r[-1][i])
               - 31824*(r[2][i] + r[-2][i])
               + 18564*(r[3][i] + r[-3][i])
               - 8568*(r[4][i] + r[-4][i])
               + 3060*(r[5][i] + r[-5][i])
               - 816*(r[6][i] + r[-6][i])
               + 153*(r[7][i] + r[-7][i])
               - 18*(r[8][i] + r[-8][i])
               + (r[9][i] + r[-9][i])
               - 48620*r0[i]
               - 1824676506132480000*r20[i]
               - 309100200138842112000*r22[i]
              ) // 6402373705728000 for i in range(len(r0))]

        r16 = [(-11440*(r[1][i] + r[-1][i])
               + 8008*(r[2][i] + r[-2][i])
               - 4368*(r[3][i] + r[-3][i])
               + 1820*(r[4][i] + r[-4][i])
               - 560*(r[5][i] + r[-5][i])
               + 120*(r[6][i] + r[-6][i])
               - 16*(r[7][i] + r[-7][i])
               + (r[8][i] + r[-8][i])
               + 12870*r0[i]
               - 4268249137152000*r18[i]
               - 527128768438272000*r20[i]
               - 51442361350668288000*r22[i]
              ) // 20922789888000 for i in range(len(r0))]

        r14 = [(3003*(r[1][i] + r[-1][i])
               - 2002*(r[2][i] + r[-2][i])
               + 1001*(r[3][i] + r[-3][i])
               - 364*(r[4][i] + r[-4][i])
               + 91*(r[5][i] + r[-5][i])
               - 14*(r[6][i] + r[-6][i])
               + (r[7][i] + r[-7][i])
               - 3432*r0[i]
               - 12204960768000*r16[i]
               - 1058170098585600*r18[i]
               - 73775500710912000*r20[i]
               - 4555411900194355200*r22[i]
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
               - 120467944396800*r18[i]
               - 5167100908569600*r20[i]
               - 208331313744153600*r22[i]
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
               - 6289809926400*r18[i]
               - 169058189664000*r20[i]
               - 4419351149875200*r22[i]
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
               - 131254905600*r18[i]
               - 2143293425280*r20[i]
               - 34682510016000*r22[i]
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
               - 771695280*r18[i]
               - 6960985920*r20[i]
               - 62711787600*r22[i]
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
               - 524280*r18[i]
               - 2097144*r20[i]
               - 8388600*r22[i]
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
               - 2*r18[i]
               - 2*r20[i]
               - 2*r22[i]
              ) // 2 for i in range(len(r0))]

        r21 = [(58786*r[1][i] 
               - 90440*r[2][i]
               + 87210*r[3][i]
               - 62016*r[4][i]
               + 33915*r[5][i]
               - 14364*r[6][i]
               + 4655*r[7][i]
               - 1120*r[8][i]
               + 189*r[9][i]
               - 20*r[10][i]
               + r[11][i]
               - 16796*r0[i]
               + 9724*r2[i]
               - 24596*r4[i]
               + 147004*r6[i]
               - 1708916*r8[i]
               + 35240284*r10[i]
               - 1237329236*r12[i]
               + 73853629564*r14[i]
               - 7850527669556*r16[i]
               + 1717351379730844*r18[i]
               - 1359124435588313876*r20[i]
               - 1272410676942417239876*r22[i]
              ) // 51090942171709440000 for i in range(len(r0))]

        r19 = [(-16796*r[1][i] 
               + 25194*r[2][i]
               - 23256*r[3][i]
               + 15504*r[4][i]
               - 7752*r[5][i]
               + 2907*r[6][i]
               - 798*r[7][i]
               + 152*r[8][i]
               - 18*r[9][i]
               + r[10][i]
               + 4862*r0[i]
               - 2860*r2[i]
               + 7436*r4[i]
               - 46420*r6[i]
               + 576236*r8[i]
               - 13098580*r10[i]
               + 532310636*r12[i]
               - 39968611540*r14[i]
               + 6350631494636*r16[i]
               - 3730771315561300*r18[i]
               - 2637991952943407764*r20[i]
               - 46833363657400320000*r21[i]
               - 734121065118879803860*r22[i]
              ) // 121645100408832000 for i in range(len(r0))]

        r17 = [(4862*r[1][i] 
               - 7072*r[2][i]
               + 6188*r[3][i]
               - 3808*r[4][i]
               + 1700*r[5][i]
               - 544*r[6][i]
               + 119*r[7][i]
               - 16*r[8][i]
               + r[9][i]
               - 1430*r0[i]
               + 858*r2[i]
               - 2310*r4[i]
               + 15258*r6[i]
               - 206790*r8[i]
               + 5386458*r10[i]
               - 272513670*r12[i]
               + 30255826458*r14[i]
               - 12765597850950*r16[i]
               - 6622557957272742*r18[i]
               - 101370917007360000*r19[i]
               - 1375210145685786630*r20[i]
               - 17172233341046784000*r21[i]
               - 201832098313986359142*r22[i]
              ) // 355687428096000 for i in range(len(r0))]

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
               - 266765571072000*r17[i]
               - 3083760849804024*r18[i]
               - 32945548027392000*r19[i]
               - 332500281299403096*r20[i]
               - 3215147584416768000*r21[i]
               - 30076927429146721464*r22[i]
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
               - 75583578470400*r17[i]
               - 643521842437836*r18[i]
               - 5269678622208000*r19[i]
               - 41890044885642492*r20[i]
               - 325386564299596800*r21[i]
               - 2481686964269990316*r22[i]
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
               - 10038995366400*r17[i]
               - 66394067988388*r18[i]
               - 430591742380800*r19[i]
               - 2750479262009668*r20[i]
               - 17360942812012800*r21[i]
               - 108550450893568228*r22[i]
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
               - 628980992640*r17[i]
               - 3275389222070*r18[i]
               - 16905818966400*r19[i]
               - 86665431465638*r20[i]
               - 441935114987520*r21[i]
               - 2244295389943190*r22[i]
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
               - 16406863200*r17[i]
               - 66398623804*r18[i]
               - 267911678160*r19[i]
               - 1078605601420*r20[i]
               - 4335313752000*r21[i]
               - 17403958407004*r22[i]
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
               - 128615880*r17[i]
               - 386371918*r18[i]
               - 1160164320*r19[i]
               - 3482590102*r20[i]
               - 10451964600*r21[i]
               - 31364282398*r22[i]
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
               - 131070*r17[i]
               - 262142*r18[i]
               - 524286*r19[i]
               - 1048574*r20[i]
               - 2097150*r21[i]
               - 4194302*r22[i]
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
              - r17[i]
              - r18[i]
              - r19[i]
              - r20[i]
              - r21[i]
              - r22[i]
               for i in range(len(r0))]
        
        return (r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,
                r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22)
        
    else:
        raise ValueError("We haven't implemented Toom-" + str(n))
        
        
        
#f= [189, 1322, 1817, 1349, 1093, 178, 1416, 646, 1031, 802, 787, 853, 33, 1557, 1619, 90]
#g= [674, 98, 846, 1892, 309, 1230, 1980, 1789, 1889, 875, 1958, 229, 1953, 146, 1842, 336]
#f = [189, 1322, 1817, 1349, 1093, 178, 1416, 646, 1031, 802, 787, 853, 33, 1557, 1619, 90] 
#g = [674, 98, 846, 1892, 309, 1230, 1980, 1789, 1889, 875, 1958, 229, 1953, 146, 1842, 336]
#f = [int(x) for x in np.random.randint(0, 100, size=40)]
#g = [int(x) for x in np.random.randint(0, 100, size=40)]
#f = [1,2,50,4]
#g = [1,6,8,7]
#print(multiply(f,g,[6]))
#print(schoolbook(f,g))


        
        
        
        
        
        
        