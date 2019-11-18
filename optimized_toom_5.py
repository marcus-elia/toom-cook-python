# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 18:36:40 2019

@author: Marcus
This is for verifying the correctness of the optimized Toom-5 formulas I
have found, before doing it in C++.
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



def optimized_toom_5(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    if len(f) > 250:
        r = {eval_list[i] : optimized_toom_5(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    else:
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    
    r0 = r[0]
        
    r8 = r['infinity']
    
    TEMP1 = [(r[1][i] + r[-1][i])//2 - r0[i] - r8[i] for i in range(len(r0))]
    
    TEMP2 = [((r[2][i] + r[-2][i])//2 
              - r0[i] 
              - 256*r8[i])//4 for i in range(len(r0))]
    
    TEMP3 = [((r[3][i] + r[-3][i])//2 - r0[i] - 6561*r8[i])//9 
             - TEMP2[i] for i in range(len(r0))]
    
    TEMP4 = [TEMP2[i] - TEMP1[i] for i in range(len(r0))]
    
    r6 = [(3*TEMP3[i] - 5*TEMP4[i])//120 for i in range(len(r0))]
    
    r4 = [(TEMP4[i] - 15*r6[i])//3 for i in range(len(r0))]
    
    r2 = [TEMP1[i] - r4[i] - r6[i] for i in range(len(r0))]
    
    TEMP5 = [(r[1][i] - r[-1][i])//2 for i in range(len(r0))]
    
    TEMP6 = [(r[2][i] - r[-2][i])//4 for i in range(len(r0))]
    
    TEMP7 = [(r[3][i] - r[-3][i])//6 for i in range(len(r0))]
    
    TEMP8 = [(r[4][i] - r0[i] - 16*r2[i] - 256*r4[i] - 4096*r6[i] - 65536*r8[i])//4
             - TEMP5[i] for i in range(len(r0))]
    
    TEMP9 = [(TEMP7[i] - TEMP5[i])//8 for i in range(len(r0))]
    
    TEMP10 = [(TEMP6[i] - TEMP5[i])//3 for i in range(len(r0))]
    
    TEMP11 = [(TEMP9[i] - TEMP10[i])//5 for i in range(len(r0))]
    
    TEMP12 = [(TEMP8[i] - 15*TEMP10[i])//180 for i in range(len(r0))]
    
    r7 = [(TEMP12[i] - TEMP11[i])//7 for i in range(len(r0))]
    
    r5 = [TEMP11[i] - 14*r7[i] for i in range(len(r0))]
    
    r3 = [TEMP10[i] - 5*r5[i] - 21*r7[i] for i in range(len(r0))]
    
    r1 = [TEMP5[i] - r3[i] - r5[i] - r7[i] for i in range(len(r0))]
    
    r_coefs = (r0, r1, r2, r3, r4, r5, r6, r7, r8)
    
    for i in range(len(r0)):
        print(i, r1[i] % 2048)
    # recombination
    k = int(np.ceil(len(f) / n))
    prod = r_coefs[0][:k]
    for j in range(1, 2*n-2):
        prod = prod + [r_coefs[j-1][k+i] + r_coefs[j][i] for i in range(k-1)]
        prod = prod + [r_coefs[j][k-1]]

    prod = prod + [r_coefs[2*n-3][k+i] + r_coefs[2*n-2][i] for i in range(k-1)]

    prod = prod + r_coefs[2*n-2][k-1:]

    return prod[:2*len(f)-1]
    



def plain_toom_5(f, g):
    n = 5
    
    fblocks = split(f, n)
    gblocks = split(g, n)
    
    # the list of evaluating numbers
    eval_list =  (0, 1, -1, 2, -2, 3, -3, 4, 'infinity')
    
    # plug the numbers in
    f_eval = evaluate_blocks_list(fblocks, eval_list)
    g_eval = evaluate_blocks_list(gblocks, eval_list)
    
    # perform the recursive multiplication
    if len(f) > 250:
        r = {eval_list[i] : plain_toom_5(f_eval[i], g_eval[i])
                for i in range(len(f_eval))}
    else:
        r = {eval_list[i] : schoolbook(f_eval[i], g_eval[i])
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

def test_optimized_toom_5():
    for poly_length in range(20, 25):
        for _ in range(100):
            f = [int(x) for x in np.random.randint(0, 2048, size=poly_length)]
            g = [int(x) for x in np.random.randint(0, 2048, size=poly_length)]
            exp = schoolbook(f, g)
            obs = optimized_toom_5(f, g)
            for i in range(len(exp)):
                if exp[i] != obs[i]:
                    print(f, g)
                    for j in range(len(exp)):
                        print(j, exp[j], obs[j], exp[j] - obs[j])
                    break
        

def compare_toom_5(min_deg, max_deg, num_trials, 
                            do_schoolbook=True, 
                            filename=None):
    """ Toggle schoolbook if you want to include schoolbook."""
    
    # the degrees (x-values)
    degrees = range(min_deg, max_deg+1)
    
    # a dictionary mapping each algorithm list to a dictionary mapping each
    # degree to the time sum for it
    optimized_deg2time_sum =  {degree:0 for degree in degrees}
    plain_deg2time_sum =  {degree:0 for degree in degrees} 
    if do_schoolbook:
        schoolbook_deg2time_sum =  {degree:0 for degree in degrees} 
                             
    # this is the outer loop so all variables will be affected equally by
    # slow outliers
    progress_time = time.time() # for the progress bar ;)
    for _ in range(num_trials):
        for degree in degrees:
            f = [int(x) for x in np.random.randint(0, 2048, degree)]
            g = [int(x) for x in np.random.randint(0, 2048, degree)]
            
            # time the optimized
            start_time = time.time()
            optimized_toom_5(f, g)
            total_time = time.time() - start_time
            optimized_deg2time_sum[degree] += total_time
            
            # time the plain
            start_time = time.time()
            plain_toom_5(f, g)
            total_time = time.time() - start_time
            plain_deg2time_sum[degree] += total_time
            
            # time schoolbook
            if do_schoolbook:
                # time the optimized
                start_time = time.time()
                schoolbook(f, g)
                total_time = time.time() - start_time
                schoolbook_deg2time_sum[degree] += total_time
                
        # progress bar for meeeeee
        if _ % (num_trials//10) == 0:
            x = _ // (num_trials//10)
            this_time = format_time(time.time()-progress_time)
            print("[" + "-"*x + " "*(10-x) + "]  " + str(this_time))
    print("[----------]")
                
    # take the average
    optimized_y_values = [optimized_deg2time_sum[degree]/num_trials for 
         degree in degrees]
    plain_y_values = [plain_deg2time_sum[degree]/num_trials for 
         degree in degrees]
    if do_schoolbook:
        schoolbook_y_values = [schoolbook_deg2time_sum[degree]/num_trials for 
         degree in degrees]
    
    # now plot
    plt.plot(degrees, optimized_y_values, lw=3, label="Optimized Toom-5")
    plt.plot(degrees, plain_y_values, lw=3, label="Plain Toom-5")
    if do_schoolbook:
        plt.plot(degrees, schoolbook_y_values, lw=3, label="Schoolbook")

    plt.xlabel("Degree", size=15)
    plt.ylabel("Average Time (seconds)", size=15)
    plt.legend(fontsize=13)
    plt.title("Running Time of Multiplication Algorithms", size=15)
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return


#compare_toom_5(735, 750, 200, do_schoolbook=False)#, filename="is_optimized_toom5_better.jpg")
f = range(1, 151)
g = [150 -x for x in range(150)]
optimized_toom_5(f,g)

