# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 10:25:54 2019

@author: Marcus
"""
import numpy as np

def schoolbook_multiply(f,g):
    """ Just plain multiplies f and g by distributing"""
    d = len(f) + len(g) - 1
    
    # create a list of zeros
    product = [0]*d
    
    for i in range(len(f)):
        for j in range(len(g)):
            product[i + j] += f[i]*g[j]
            
    return product

def poly_add(f,g):
    return [f[i] + g[i] for i in range(len(f))]

def poly_subtract(f,g):
    return [f[i] - g[i] for i in range(len(f))]


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

def single_karatsuba(f,g):
    """ Performs a single iteration of Karatsuba multiplication and does
        schoolbook for the three smaller mults."""
    (f0, f1) = split(f,2)
    (g0, g1) = split(g,2)
    
    r0 = schoolbook_multiply(f0, g0)
    r2 = schoolbook_multiply(f1, g1)
    temp_prod = schoolbook_multiply(poly_add(f0,f1), poly_add(g0,g1))
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1]   

    
def karatsuba_multiply(f,g, minsize=10):
    """ Performs recursive karatsuba multiplication until minsize"""
    if len(f) < minsize:
        return schoolbook_multiply(f,g)
    
    (f0, f1) = split(f,2)
    (g0, g1) = split(g,2)
    
    r0 = karatsuba_multiply(f0, g0, minsize=minsize)
    r2 = karatsuba_multiply(f1, g1, minsize=minsize)
    temp_prod = karatsuba_multiply(poly_add(f0,f1), poly_add(g0,g1), minsize=minsize)
    r1 = poly_subtract(poly_subtract(temp_prod, r2), r0)
    
    k = int(np.ceil(len(f) / 2))
    product = r0[:k]
    product = product + [r0[k+i] + r1[i] for i in range(k-1)]
    product = product + [r1[k-1]]
    product = product + [r1[k+i] + r2[i] for i in range(k-1)]
    product = product + r2[k-1:]
    return product[:2*len(f)-1] 
    
    # this is the most confusing thing ever
    # deprecated
    #if len(f) % 2 == 0:
    #    k = len(f)//2
    #    part1 = r0[:k] 
    #    part2 = [r0[k+i] + r1[i] for i in range(k-1)] 
    #    part3 = [r1[k-1]]
    #    part4 = [r1[k+i] + r2[i] for i in range(k-1)]
    #    part5 = r2[k-1:]
    #    return part1 + part2 + part3 + part4 + part5
    #else:
    #    k = len(f)//2
    #    part1 = r0[:k+1]
    #    part2 = [r0[k+1+i] + r1[i] for i in range(k)]
    #    part3 = [r1[k]]
    #    part4 = [r1[k+i+1] + r2[i] for i in range(0,k)]
    #    part5 = r2[k:4*k+1 - 2*k-2]
    #    return part1 + part2 + part3 + part4 + part5
    
#f = [8,7,6,5,4,3,2]
#g = [1,2,3,4,5,6,7]
#print(schoolbook_multiply(f,g))
#print(single_karatsuba(f,g))