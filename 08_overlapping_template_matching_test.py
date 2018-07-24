
from __future__ import print_function

import math
import scipy.special as ss

def lgamma(x):
    return math.log(gamma(x))
    
def Pr(u, eta):
    if ( u == 0 ):
        p = math.exp(-eta)
    else:
        sum = 0.0
        for l in range(1,u+1):
            sum += math.exp(-eta-u*math.log(2)+l*math.log(eta)-lgamma(l+1)+lgamma(u)-lgamma(l)-lgamma(u-l+1))
        p = sum
    return p

def test(input, n, blen=6):
    
    m = 10
    # Build the template B as a random list of input
    B = [1 for x in range(m)]
    
    N = 968
    K = 5
    M = 1062
    
    if len(input) < (M*N):
        print("Insufficient data. %d bit provided. 1,028,016 input required" % len(input))
        return 0.0, 0.0, False
    
    blocks = list() # Split into N blocks of M input
    for i in range(N):
        blocks.append(input[i*M:(i+1)*M])

    # Count the distribution of matches of the template across blocks: Vj
    v=[0 for x in range(K+1)] 
    for block in blocks:
        count = 0
        for position in range(M-m):
            if block[position:position+m] == B:
                count += 1
            
        if count >= (K):
            v[K] += 1
        else:
            v[count] += 1

    chisq = 0.0  # Compute Chi-Square

    pi = [0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865] # From STS
    piqty = [int(x*N) for x in pi]
    
    lambd = (M-m+1.0)/(2.0**m)
    eta = lambd/2.0
    sum = 0.0
    for i in range(K): #  Compute Probabilities
        pi[i] = Pr(i, eta)
        sum += pi[i]

    pi[K] = 1 - sum;

    sum = 0    
    chisq = 0.0
    for i in range(K+1):
        chisq += ((v[i] - (N*pi[i]))**2)/(N*pi[i])
        sum += v[i]
        
    p = ss.gammaincc(5.0/2.0, chisq/2.0) # Compute P value

    print("  B = ",B)
    print("  m = ",m)
    print("  M = ",M)
    print("  N = ",N)
    print("  K = ",K)
    print("  model = ",piqty)
    print("  v[j] =  ",v) 
    print("  chisq = ",chisq)
    
    success = ( p >= 0.01)
    return [chisq, p, success]
