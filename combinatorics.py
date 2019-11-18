# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:09:30 2019

@author: Marcus

This file is for helping me to prove the identity
"""

import itertools

def standardize(string, permutation):
    """ string is a list of positive integers, and so is permutation.
        permutation is a permutation of 1..n. probably either 1,2,...,n
        or n, n-1,...,1. string is a guide of how to do the permutation.
        So 1,2,1,2 according to 4,3,2,1 = 4,2,3,1"""
    if len(string) != len(permutation):
        raise ValueError("They have to be same size, not {} and {}".format(len(string), len(permutation)))
        
    result = [0 for _ in range(len(string))]
    permutation_index = 0
    min_value = 1
    while permutation_index < len(permutation):
        for i in range(len(string)):
            current_int = string[i]
            if current_int == min_value:
                result[i] = permutation[permutation_index]
                permutation_index += 1
        min_value += 1
    return tuple(result)

def standardize_set(m):
    results = []
    for string in itertools.product(range(1, m+1), repeat=2*m):
        results.append(standardize(string, range(1, 2*m + 1)))
    return results

def standardize_dict(m):
    str2per = {}
    for string in itertools.product(range(1, m+1), repeat=2*m):
        str2per[string] = standardize(string, range(1, 2*m + 1))
    return str2per

def reverse_list(L):
    n = len(L)
    return [L[n - i - 1] for i in range(n)]

def reverse_set(m):
    return [reverse_list(L) for L in standardize_set(m)]


def check_reverses(m):
    L = standardize_set(m)
    R = reverse_set(m)
    for permutation in itertools.permutations(range(1,2*m+1)):
        if not tuple(permutation) in L and not tuple(reverse_list(permutation)) in L:
            print(permutation)
            return False
    return True
        
def num_descents(permutation):
    count = 0
    for i in range(len(permutation) - 1):
        if permutation[i] > permutation[i + 1]:
            count += 1
    return count

# not useful
def descent_locations(string):
    """ Returns a list of the indices where a descent happens."""
    locs = []
    for i in range(len(string) - 1):
        if string[i] > string[i + 1]:
            locs.append(i)
    return locs

def des_asc_locations(string):
    """ Puts a 1 at every ascent location, -1 at every descent location,
        and 0 elsewhere."""
    locs = [0 for _ in range(len(string) - 1)]
    for i in range(len(string) - 1):
        if string[i] > string[i + 1]:
            locs[i] = -1
        elif string[i] < string[i+1]:
            locs[i] = 1
    return locs

def num_inv_descents(permutation):
    count = 0
    for i in range(1, len(permutation)):
        for j in range(0, i):
            if permutation[j] == permutation[i] + 1:
                count += 1
    return count



def get_map_results(m):
    results = []
    permutation = list(range(1, 2*m + 1))
    for string in itertools.product(range(1, m+1), repeat=2*m):
        results.append(standardize(string, permutation))
    #permutation.reverse()
    #for string in itertools.product(range(1, m+1), repeat=2*m):
    #    results.append(standardize(string, permutation))
    return results

def get_map_dictionary(m):
    results = {}
    permutation = list(range(1, 2*m + 1))
    for string in itertools.product(range(1, m+1), repeat=2*m):
        results[string] = [standardize(string, permutation)]
    #permutation.reverse()
    #for string in itertools.product(range(1, m+1), repeat=2*m):
    #    results[string].append(standardize(string, permutation))
    return results

def get_inverse_dictionary(m):
    map_dictionary = get_map_dictionary(m)
    inverse_dictionary = {}
    for string in map_dictionary:
        for permutation in map_dictionary[string]:
            try:
                inverse_dictionary[permutation].append(string)
            except KeyError:
                inverse_dictionary[permutation] = [string]
    return inverse_dictionary

def print_repeats(m):
    inverse_dictionary = get_inverse_dictionary(m)
    for permutation in inverse_dictionary:
        if len(inverse_dictionary[permutation]) > 1:
            print(permutation, inverse_dictionary[permutation])
            print("\n")
    return

def get_repeat_dictionary(m):
    """ Makes a dictionary that maps a number of repetitions to the number of
        permutations that are repeated that many times"""
    inverse_dictionary = get_inverse_dictionary(m)
    repeat_dictionary = {}
    for permutation in inverse_dictionary:
        num_repeats = len(inverse_dictionary[permutation])
        try:
            repeat_dictionary[num_repeats] += 1
        except KeyError:
            repeat_dictionary[num_repeats] = 1
    return repeat_dictionary

def check_surjection(m):
    """ Checks that every permutation appears in our list of
        permutations resulting from standardizing strings."""
    outputs = set(get_map_results(m))
    desired_range = set(itertools.permutations(range(1, 2*m + 1)))
    return outputs == desired_range
    
def is_ascending(L):
    """ Checks if list L is ascending"""
    for i in range(len(L) - 1):
        if L[i+1] < L[i]:
            return False
    return True

def num_ascents(permutation):
    count = 0
    for i in range(len(permutation) - 1):
        if permutation[i] < permutation[i + 1]:
            count += 1
    return count


def split_from_right(chain, character):
    """ chain must be ascending, and character must be >= first elmt of chain
        and < last elmt of chain. So
        split_from_right([1,2,3,4],1) = [[1], [2,3,4]]
        split_from_right([1,2,3,4],3) = [[1,2,3], [4]]
        split_from_right([1,2,3,4],4) = error
        split_from_right([1,3],2) = [[1],[3]]
        split_from_right([1,1,2],1) = [[1,1], [2]]"""
    if character >= chain[-1]:
        raise ValueError("can't right split {} on {}".format(chain, character))
    i = len(chain) - 1
    while i >= 0:
        if chain[i] <= character:
            return [chain[:i+1], chain[i+1:]]
        i -= 1
    # in case the character is strictly less than the last element, but
    # wasn't <= any previous element
    if chain[-1] > character:
        return [chain[:-1], [chain[-1]]]
    raise ValueError("can't right split {} on {}".format(chain, character))
    return
    
def split_from_left(chain, character):
    """ chain must be ascending, and character must be > first elmt of chain
        and <= last elmt of chain. So
        split_from_left([1,2,3,4],1) = error
        split_from_left([1,2,3,4],3) = [[1,2], [3,4]]
        split_from_left([1,2,3,4],4) = [[1,2,3],[4]]
        split_from_left([1,3], 2) = [[1], [3]]
        split_from_left([1,1,2],1) = error"""
    if character <= chain[0]:
        raise ValueError("can't left split {} on {}".format(chain, character))
    i = 0
    while i < len(chain) - 1:
        if chain[i] < character <= chain[i+1]:
            return [chain[:i+1], chain[i+1:]]
        i += 1
    # in case the character is strictly greater than the first element, but
    # wasn't >= any other element
    if chain[0] < character:
        return[[chain[0]], chain[1:]]
    raise ValueError("can't left split {} on {}".format(chain, character))
    return



def split_one_step(chain_list):
    """ This attempts to perform one step of splitting, and returns the result
        If no splitting needs to be done it returns the original chain list."""
    # reconstruct the original string
    string = []
    for chain in chain_list:
        string = string + chain
        
    # do the first chain
    first_chain = chain_list[0]
    min_a = first_chain[0]
    max_a = first_chain[-1]
    for chain in chain_list[1:]:
        for character in chain:
            if character >= min_a and character < max_a:
                return split_from_right(first_chain, character) + chain_list[1:]
            
    # do all the middle chains
    for i in range(1, len(chain_list) - 1):
        current_chain = chain_list[i]
        min_a = current_chain[0]
        max_a = current_chain[-1]
        
        # the left side
        for chain in chain_list[:i]:
            for character in chain:
                if character > min_a and character <= max_a:
                    return chain_list[:i] + split_from_left(current_chain, character) + chain_list[i+1:]
        
        # the right side
        for chain in chain_list[i+1:]:
            for character in chain:
                if character >= min_a and character < max_a:
                    return chain_list[:i] + split_from_right(current_chain, character) + chain_list[i+1:]
    
    # do the last chain
    last_chain = chain_list[-1]
    min_a = last_chain[0]
    max_a = last_chain[-1]
    for chain in chain_list[:-1]:
        for character in chain:
            if character > min_a and character <= max_a:
                return chain_list[:-1] + split_from_left(last_chain, character)
            
    # if no split was found
    return chain_list

    
def ascpr(string):
    """ Returns the (improved) ascending chain priority of a list.
        For example, ascpr([1,3,1,2,2,4]) = [(1,1), (3,1), (2,3), (4,1)]"""
    # special case
    if is_ascending(string):
        return[(1,len(string))]
    
    chain_list = [] 

    current_chain = [string[0]]       
    for i in range(len(string) - 1):
        
        # if the next character is in the chain
        if string[i] <= string[i+1]:
            
            # add the character to the chain
            current_chain.append(string[i+1])
            
        # otherwise, the chain ends
        else:
            chain_list.append(current_chain)
            current_chain = [string[i+1]]
            
        
    chain_list.append(current_chain)
    
    # now we have to go back and split things up
    X = split_one_step(chain_list)
    while X != chain_list:
        chain_list = X
        X = split_one_step(chain_list)
    
    # now we find the priorities
    index2priority = {} # map the index of the chain to its priority
    min_element = 1
    priority = 1
    while len(index2priority) < len(chain_list):
        for i in range(len(chain_list)):
            chain = chain_list[i]
            if chain[0] == min_element:
                index2priority[i] = priority
                priority += 1
        min_element += 1
    priority_list = []
    for i in range(len(chain_list)):
        chain = chain_list[i]
        priority_list.append((index2priority[i], len(chain)))
    return priority_list

def print_repeats_ascpr(m):
    inverse_dictionary = get_inverse_dictionary(m)
    for permutation in inverse_dictionary:
        print(len(inverse_dictionary[permutation]), 
              ascpr(inverse_dictionary[permutation][0]),
              ascpr_inv_descents(inverse_dictionary[permutation][0]))
        break
        #print("\n")
    return

def ascpr_inv_descents(string):
    """ Returns the number of inverse descents in an ascpr."""
    asc = ascpr(string)
    return num_inv_descents([a[0] for a in asc])

def count_permutations_with_num_inv_descents(m):
    """ Returns a dictionary mapping the number of inv descents a permutation
        has to the number of permutations with that many"""
    num2count = {}
    for permutation in itertools.permutations(range(2*m)):
        d = num_inv_descents(permutation)
        try:
            num2count[d] += 1
        except KeyError:
            num2count[d] = 1
    return num2count

def num_ascents(permutation):
    count = 0
    for i in range(len(permutation) - 1):
        if permutation[i] < permutation[i+1]:
            count += 1
    return count

def get_permutations_with_target_ascents(m, target):
    permutations = []
    for permutation in itertools.permutations(range(1,2*m+1)):
        if num_ascents(permutation) == target:
            permutations.append(permutation)
    return permutations

def get_permutations_with_target_inv_descents(m, target):
    permutations = []
    for permutation in itertools.permutations(range(1,2*m+1)):
        if num_inv_descents(permutation) == target:
            permutations.append(permutation)
    return permutations

def inverse_standardization(sigma):
    """ Returns the minimal string that maps to sigma"""
    s = [0]*len(sigma)
    character = 1
    previous_location = 0
    for i in range(1, len(sigma) + 1):
        current_location = sigma.index(i)
        
        # if we have an inverse descent
        if current_location < previous_location:
            character += 1
        s[current_location] = character
        previous_location = current_location
    return s

def increment(s, i):
    """ Returns the result of incrementing s at index i. So
        1,2,3,1,2,3,4,3 incremented at 2 becomes 1,2,4,1,2,4,5,4"""
    incremented = list(s)
    turning_point = s[i]
    for j in range(0, i):
        if incremented[j] > turning_point:
            incremented[j] += 1
    for j in range(i, len(incremented)):
        if incremented[j] >= turning_point:
            incremented[j] += 1
    return incremented

def increment_list(s, indices):
    incremented = list(s)
    for i in indices:
        incremented = increment(incremented, i)
    return incremented

def increment_list_mod(s, indices, m):
    """ Does it somewhat mod m, except we use m and not 0"""
    incremented = list(s)
    for i in indices:
        incremented = increment(incremented, i)
    return [x % m if x % m != 0 else m for x in incremented]

def get_equivalence_class(s):
    """ Takes the minimal string, s, and returns the list of strings
        that can be made by adding a single increment."""
    repeated_strings = []
    for i in range(0, len(s)):
        repeated_strings.append(increment(s, i))
    return repeated_strings

def get_double_subtracted_strings(m):
    string2subtracts = {}  # maps each string to the number of times
                            # it gets subtracted off
    
    for s in itertools.product(range(1,m), repeat=2*m):
        for string in get_equivalence_class(s):
            try:
                string2subtracts[tuple(string)] += 1
            except KeyError:
                string2subtracts[tuple(string)] = 1
    #return [(string, string2subtracts[string]) for string in
    #        string2subtracts if string2subtracts[string] > 1]
    return string2subtracts
    
def count_subtracted_strings_dictionary(m):
    """ Returns a dictionary mapping the number of times a string gets
        subtracted off in the first exclusion to the number of strings
        that are subtracted that many times. So for m=2, we expect to get"""
    string2subtracts = get_double_subtracted_strings(m)
    subtracts2numstrings = {}
    for string in string2subtracts:
        n = string2subtracts[string]
        try:
            subtracts2numstrings[n] += 1
        except KeyError:
            subtracts2numstrings[n] = 1
    return subtracts2numstrings

def get_incremented_strings(string, m, j):
    incremented_strings = []
    if j == m-1:
        choose_amounts = (0, 2*m-3, 2*m-2)
    elif j > 1:
        choose_amounts = (m-j-2, m-j-1, m+j-1, m-j)
    
    for choose_amount in choose_amounts:
        for indices in itertools.combinations(range(2*m - 2), choose_amount):
            incremented_strings.append(tuple(increment_list_mod(string, indices, m)))
    return incremented_strings

def get_set_of_strings(m, j):
    strings = []
    for string in itertools.product(range(1, j+1), repeat=2*m-2):
        strings = strings + get_incremented_strings(string, m, j)
    return strings




#for p in standardize_set(4):
#    if num_inv_descents(p) > 3:
#        print('bad')
#print()
#for p in reverse_set(4):
#    if num_inv_descents(p) < 4:
#        print('bad')
#print(check_reverses(4))

# ===================================
#
#  Trying to prove the 0 identities
#
# ===================================

def toInt(coefs):
    """ Converts the list [1,2,3] to the integer 123"""
    i = 1
    integer = 0
    multiplier = 1
    while i < len(coefs) + 1:
        integer += multiplier * coefs[-i]
        i += 1
        multiplier *= 10
    return integer

def getIntegers(m, j):
    """ Returns the integers you get when you take strings of length 2
        from {1,2, ..., j}, make them integers, and multiply them by
        the integers you get from the combinations of 2m choose m-j"""
    intList = []
    for i in itertools.product(range(1, j+1), repeat=2):
       for k in itertools.combinations(range(1,2*m+1), m-j):
           if toInt(k) == 0:
               intList.append(toInt(i))
           else:
               intList.append(toInt(i)*toInt(k))
    return intList

print(getIntegers(3, 3))


