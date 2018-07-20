from __future__ import print_function

import math
from gamma_functions import *

def int2patt(n,m):
    pattern = list()
    for i in range(m):
        pattern.append((n >> i) & 1)
    return pattern
    
def countpattern(patt,input,n):
    thecount = 0
    for i in range(n):
        match = True
        for j in range(len(patt)):
            if patt[j] != input[i+j]:
                match = False
        if match:
            thecount += 1
    return thecount

def psi_sq_mv1(m, n, padded_input):
    counts = [0 for i in range(2**m)] 
    for i in range(2**m):
        pattern = int2patt(i,m)
        count = countpattern(pattern,padded_input,n)
        counts.append(count)
        
    psi_sq_m = 0.0
    for count in counts: 
        psi_sq_m += (count**2)
    psi_sq_m = psi_sq_m * (2**m)/n 
    psi_sq_m -= n
    return psi_sq_m            
         
def serial_test(input,patternlen=None):
    

    n = len(input) - 1

    if n%2 != 0:
        print(n)
        exit()

    if patternlen != None:
        m = patternlen  
    else:  
        m = int(math.floor(math.log(n,2)))-2
    
        if m < 4:
            print("Error. Not enough data for m to be 4")
            return False,0,None
        m = 4

    # Step 1
    padded_input=input+input[0:m-1]

    padded_input = int('1001',2)

    print("padded input = " +str(padded_input))
    
    # Step 2
    psi_sq_m   = psi_sq_mv1(m, n, padded_input)
    psi_sq_mm1 = psi_sq_mv1(m-1, n, padded_input)
    psi_sq_mm2 = psi_sq_mv1(m-2, n, padded_input)    
    
    delta1 = psi_sq_m - psi_sq_mm1
    delta2 = psi_sq_m - (2*psi_sq_mm1) + psi_sq_mm2
    
    p1 = gammaincc(2**(m-2),delta1/2.0)
    p2 = gammaincc(2**(m-3),delta2/2.0)
        
    # print("  psi_sq_m   = ",psi_sq_m)
    # print("  psi_sq_mm1 = ",psi_sq_mm1)
    # print("  psi_sq_mm2 = ",psi_sq_mm2)
    # print("  delta1     = ",delta1)
    # print("  delta2     = ",delta2)  
    # print("  p1         = ",p1)
    # print("  p2         = ",p2)
     
    success = (p1 >= 0.01) and (p2 >= 0.01)
    return [psi_sq_m, psi_sq_mm1, psi_sq_mm2, delta1, delta2, p1, p2, success]
    
