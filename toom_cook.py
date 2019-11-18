# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 10:45:48 2018

@author: Marcus

In this document, I am primarily learning how to implement
the Toom Cook algorithm for multiplying large polynomials.
"""
import numpy as np


# ==========================================
#
# TOOM COOK
#
# ==========================================

def schoolbook_multiply(f,g):
    """ Just plain multiplies f and g by distributing"""
    d = len(f) + len(g) - 1
    
    # create a list of zeros
    product = [0]*d
    
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] += f[i]*g[j]
            
    return product

def remove_excess_zeros(f):
    while len(f) > 1 and f[-1] == 0:
        f = f[:-1]
    return f

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

def multiply_corresponding(blocks1, blocks2, n=3,  minsize=32):
    """ the blocks are lists of lists, each list being the coefficients
        of a polynomial. This uses multiplication to
        multiply the corresponding entries."""
    answer = []
    for i in range(len(blocks1)):
        answer.append(toom_cook_multiply(blocks1[i], blocks2[i], n=n, minsize=minsize))
    return answer

# =======================================
# ========= Fraction Stuff ==============
# =======================================
class Fraction():
    
    def __init__(self, num, den):
        # reduce this garbage
        while num % 2 == 0 and den % 2 == 0:
            num = num // 2
            den = den // 2
        while num % 3 == 0 and den % 3 == 0:
            num = num // 3
            den = den // 3
        while num % 5 == 0 and den % 5 == 0:
            num = num // 5
            den = den // 5
        self.numer = num
        self.denom = den
    
    def __repr__(self):
        if self.denom == 1:
            return str(self.numer)
        return str(self.numer) + "/" + str(self.denom)
        
    def to_int(self):
        """ Returns an int equal to the fraction. Throws an error if this
            is not possible."""
        if self.numer % self.denom == 0:
            return self.numer//self.denom
        else:
            raise ValueError("{}/{} is not equal to an integer".format(self.numer, self.denom))
    
    def scalar_mult(self, scalar):
        """ Returns a new fraction that is this fraction
            multiplied by the scalar."""
        return Fraction(self.numer * scalar, self.denom * scalar)
    
    
def add_fractions(a, b):
    """ Adds a and b as fractions. Converts them to fractions
        if they are ints."""
    if type(a) == int:
        a = Fraction(a, 1)
    if type(b) == int:
        b = Fraction(b, 1)
    
    den = a.denom * b.denom
    num = b.denom * a.numer + a.denom * b.numer
    
    return Fraction(num, den)
    
def sub_fractions(a, b):
    """ Subtracts b from a as fractions. Converts them to fractions
        if they are ints."""
    if type(a) == int:
        a = Fraction(a, 1)
    if type(b) == int:
        b = Fraction(b, 1)
        
    den = a.denom * b.denom
    num = b.denom * a.numer - a.denom * b.numer
    return Fraction(num, den)

def mult_fractions(a, b):
    """ Multiplies a and b as fractions. Converts them to
        fractions if they are ints."""
    if type(a) == int:
        a = Fraction(a, 1)
    if type(b) == int:
        b = Fraction(b, 1)
    
    return Fraction(a.numer*b.numer, a.denom*b.denom)

# =======================================
# ======== End Fraction Stuff ===========
# =======================================

def matrix_times_vector(mat, vec):
    """ Multiplies mat by vec. Assumes some entries of mat are fractions, and
        the rest are ints. But each entry of vec is a list of ints,
        representing a polynomial. """
    if len(mat[0]) != len(vec):
        raise ValueError("Can't multiply this matrix and vector")
        
    # the number of coefficients in each polynomial
    n = len(vec[0])
    
    # initialize each entry of the resulting vector to be
    # n zeros, but each zero is a fraction.
    answer = [[Fraction(0,1) for _ in range(n)] for _ in range(len(vec))]
    
    for j in range(len(mat)):
        for i in range(len(vec)):
            for coef in range(n):
                mat_x_vec = mult_fractions(mat[j][i], vec[i][coef])
               
                answer[j][coef] = add_fractions(answer[j][coef], mat_x_vec)
   
    # convert the fraction objects to ints
    for i in range(len(answer)):
        for j in range(n):
            answer[i][j] = answer[i][j].to_int()
    return answer



def optimized_toom_3(f, g, minsize=20):
    """ An optimized version that only applies to Toom-3. Avoids explicitly
        using a matrix."""
    if len(f) <= minsize:
        return schoolbook_multiply(f,g)
    
     # split the polynomials into blocks
    fblocks = split(f, 3)
    gblocks = split(g, 3)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    r = [optimized_toom_3(f_eval[i], g_eval[i], minsize=minsize) for i in range(len(f_eval))]

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


def make_toom_matrix(n=3):
    """ Makes the matrix to be used in Toom-Cook"""
    if n == 2:
        row0 = [1,0,0]
        row1 = [-1,1,-1]
        row2 = [0,0,1]
        
        return [row0, row1, row2]
    
    if n == 3:
        row0 = [1,0,0,0,0]
        row1 = [Fraction(1,2), Fraction(1,3), -1, Fraction(1,6), -2]
        row2 = [-1, Fraction(1,2), Fraction(1,2), 0, -1]
        row3 = [Fraction(-1,2), Fraction(1,6), Fraction(1,2), Fraction(-1,6), 2]
        row4 = [0,0,0,0,1]

        return [row0, row1, row2, row3, row4]
    
    if n == 4:
        row0 = [1,0,0,0,0,0,0]
        row1 = [Fraction(-1,3), 1, Fraction(-1,2), Fraction(-1,4), Fraction(1,20), Fraction(1,30), -12]
        row2 = [Fraction(-5,4), Fraction(2,3), Fraction(2,3), Fraction(-1,24), Fraction(-1,24), 0, 4]
        row3 = [Fraction(5,12), Fraction(-7,12), Fraction(-1,24), Fraction(7,24), Fraction(-1,24), Fraction(-1,24), 15]
        row4 = [Fraction(1,4), Fraction(-1,6), Fraction(-1,6), Fraction(1,24), Fraction(1,24), 0, -5]
        row5 = [Fraction(-1,12), Fraction(1,12), Fraction(1,24), Fraction(-1,24), Fraction(-1,120), Fraction(1,120), -3]
        row6 = [0,0,0,0,0,0,1]
        
        return [row0, row1, row2, row3, row4, row5, row6]
    
def optimized_toom_4(f, g, minsize=20):
    """ An optimized version that only applies to Toom-4. Avoids explicitly
        using a matrix."""
    if len(f) <= minsize:
        return schoolbook_multiply(f,g)
    
     # split the polynomials into blocks
    fblocks = split(f, 4)
    gblocks = split(g, 4)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i]:optimized_toom_4(f_eval[i], g_eval[i], minsize=minsize) for i in range(len(f_eval))}
    
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
    
def make_toom_evaluation_list(n=3):
    """ Returns the list of x-values that we plug in to the 
        polynomials for toom-cook. Use 'infinity' for infinity"""
    if n == 2:
        return [0, 1, 'infinity']
    
    if n == 3:
        return [0, 1, -1 ,-2, 'infinity']
    
    if n == 4:
        return [0, 1, -1, 2, -2, 3, 'infinity']
    
def toom_cook_multiply(f, g, n=3, minsize=32):
    """ Performs toom-cook multiplication. If the degree is small enough,
        just do schoolbook multiplication."""    
    
    if len(f) <= minsize:
        return schoolbook_multiply(f,g)
    
   
    if len(f) != len(g):
        raise ValueError("I'm not smart enough to handle polynomials of different sizes.")
    
    # split the polynomials into blocks
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_toom_evaluation_list(n=n)
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this is where the recursion happens: multiply smaller polynomials
    vector = multiply_corresponding(f_eval, g_eval, n=n, minsize=minsize)
    
    # this just constructs a predetermined inverse matrix to multiply
    matrix = make_toom_matrix(n=n)
    
    resultant_vector = matrix_times_vector(matrix, vector)
    # combine the coefficient polynomials into the product
    product = [0 for _ in range(len(f) + len(g) - 1)]
    power_of_x = 0
    for coefficient_polynomial in resultant_vector:
        for i in range(len(coefficient_polynomial)):
            
            # need this if statement to avoid index errors
            # alternatively, make a special case for the last couple
            # of polynomial coefficients
            if i + power_of_x < len(product):
                product[i + power_of_x] += coefficient_polynomial[i]
        power_of_x += int(np.ceil(len(f) / n))
           
    return product 

def single_toom_3(f,g):
    """ Performs a single (optimized) iteration of Toom-3 """
     # split the polynomials into blocks
    fblocks = split(f, 3)
    gblocks = split(g, 3)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    r = [schoolbook_multiply(f_eval[i], g_eval[i]) for i in range(len(f_eval))]

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

def single_toom_4(f,g):
    """ Performs a single (optimized) iteration of Toom-4 """
    # split the polynomials into blocks
    fblocks = split(f, 4)
    gblocks = split(g, 4)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i]) for i in range(len(f_eval))}
    
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

def single_toom_5(f,g):
    """ Performs a single (optimized) iteration of Toom-5 """
    # split the polynomials into blocks
    fblocks = split(f, 5)
    gblocks = split(g, 5)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # this contains r(0), r(1), r(-1), r(2), r(infinity)
    r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i])
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


def optimized_toom_5(f, g, minsize=20):
    """ An optimized version that only applies to Toom-5. Avoids explicitly
        using a matrix."""
    if len(f) <= minsize:
        return schoolbook_multiply(f,g)
    
     # split the polynomials into blocks
    fblocks = split(f, 5)
    gblocks = split(g, 5)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i]:optimized_toom_5(f_eval[i], g_eval[i], minsize=minsize)
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


def single_toom_6(f,g):
    """ Performs a single (optimized) iteration of Toom-6 """
    # split the polynomials into blocks
    fblocks = split(f, 6)
    gblocks = split(g, 6)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    
    r = {eval_list[i]:schoolbook_multiply(f_eval[i], g_eval[i])
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



def optimized_toom_6(f, g, minsize=20):
    """ An optimized version that only applies to Toom-6. Avoids explicitly
        using a matrix."""
    if len(f) <= minsize:
        return schoolbook_multiply(f,g)
    
     # split the polynomials into blocks
    fblocks = split(f, 6)
    gblocks = split(g, 6)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, 'infinity')
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    r = {eval_list[i]:optimized_toom_6(f_eval[i], g_eval[i], minsize=minsize)
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

#f = [1,2,3,4,5,6,7,8,9,10,11,12]
#g = [12,11,10,9,8,7,6,5,4,3,2,1]
#print(schoolbook_multiply(f,g))
#print(optimized_toom_6(f,g, minsize=5))
#print(schoolbook_multiply(f,g))
#print(optimized_toom_4(f, g, minsize=5))
#print(toom_cook_multiply(f,g,n=4,minsize=5))


