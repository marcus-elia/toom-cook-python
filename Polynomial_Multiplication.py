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
    product = [0 for _ in range(d)]
    
    for i in range(len(f)):
        for j in range(len(g)):
            product[(i + j) % d] += f[i]*g[j]
            
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
    while len(f) % num_blocks != 0:
        f.append(0)
    block_length = len(f) // num_blocks
    index = 0
    while index + block_length < len(f):
        blocks.append(f[index:index+block_length])
        index += block_length
    blocks.append(f[index:])
    return blocks

def evaluate_blocks(blocks, value):
    """ blocks is a list of lists, each list is the coefficients of a
        polynomial. But each list a coefficient. For example, if blocks is
        [[1,2],[3,4],[5,6]] and value is -2, we return
        [1,2] + [-6,-8] + [20,24] = [15, 18].  If the value is infinity,
        we return the leading coefficient."""
        
    if value == 'infinity':
        return blocks[-1]
    
    # initialize an empty of the right length
    answer = [0 for _ in range(len(blocks[0]))]
    
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

def multiply_corresponding(blocks1, blocks2):
    """ the blocks are lists of lists, each list being the coefficients
        of a polynomial. This uses schoolbook multiplication to
        multiply the corresponding entries."""
    answer = []
    for i in range(len(blocks1)):
        answer.append(schoolbook_multiply(blocks1[i], blocks2[i]))
    return answer

class Fraction():
    
    def __init__(self, num, den):
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

def make_toom_matrix(n=3):
    """ Makes the matrix to be used in Toom-Cook"""
    if n == 3:
        row0 = [1,0,0,0,0]
        row1 = [Fraction(1,2), Fraction(1,3), -1, Fraction(1,6), -2]
        row2 = [-1, Fraction(1,2), Fraction(1,2), 0, -1]
        row3 = [Fraction(-1,2), Fraction(1,6), Fraction(1,2), Fraction(-1,6), 2]
        row4 = [0,0,0,0,1]

        return [row0, row1, row2, row3, row4]
    
def make_toom_evaluation_list(n=3):
    """ Returns the list of x-values that we plug in to the 
        polynomials for toom-cook. Use 'infinity' for infinity"""
    if n == 3:
        return [0, 1, -1 ,-2, 'infinity']
    
def toom_cook_multiply(f,g,n=3):
    """ Performs toom-cook multiplication"""
    
    # split the polynomials into blocks
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  make_toom_evaluation_list(n=n)
    
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    vector = multiply_corresponding(f_eval, g_eval)
    matrix = make_toom_matrix(n=n)
    
    resultant_vector = matrix_times_vector(matrix, vector)
    
    # combine the coefficient polynomials into the product
    product = [0 for _ in range(len(f) + len(g) - 1)]
    power_of_x = 0
    for coefficient_polynomial in resultant_vector:
        for i in range(len(coefficient_polynomial)):
            product[i + power_of_x] += coefficient_polynomial[i]
        power_of_x += int(np.ceil(len(f) / n))
           
    return remove_excess_zeros(product)   


for _ in range(10):
    f = [int(x) for x in np.random.randint(0, 101, size=500)]
    g = [int(x) for x in np.random.randint(0, 101, size=500)]
    true_prod = schoolbook_multiply(f,g)
    prod = toom_cook_multiply(f,g)
    if true_prod != prod:
        print(true_prod)
        print(prod)



print(toom_cook_multiply([1,2,3,4,5,6], [8,7,6,5,4,3]))