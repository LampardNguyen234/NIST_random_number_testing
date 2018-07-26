import math
import scipy.special as ss

def berelekamp_massey(input):
    n = len(input)
    b = [0 for x in input]  #initialize b and c arrays
    c = [0 for x in input]
    b[0] = 1
    c[0] = 1
    
    L = 0
    m = -1
    N = 0
    while (N < n):
        #compute discrepancy
        d = int(input[N],2)
        for i in range(1,L+1):
            if input[N-i] == '1':
                t = 1
            else:
                t = 0
            d = d ^ (c[i] & t)
            # print("d = " + str(d))
        if (d != 0):  # If d is not zero, adjust poly
            t = c[:]
            for i in range(0,n-N+m):
                c[N-m+i] = c[N-m+i] ^ b[i] 
            if (L <= (N/2)):
                L = N + 1 - L
                m = N
                b = t 
        N = N +1
    # Return length of generator and the polynomial
    # print("Done")
    return L , c[0:L]
    
def test(input, n, patternlen=None):
    n = len(input)
    # Step 1. Choose the block size
    if patternlen != None:
        M = patternlen  
    else: 
        if n < 1000000:
            print("Error. Need at least 10^6 input")
            return [0]*9
        M = 512
    K = 6 
    N = int(math.floor(n/M))
    
    # Step 2 Compute the linear complexity of the blocks
    LC = list()
    for i in range(N):
        # print("i = " + str(i))
        x = input[(i*M):((i+1)*M)]
        LC.append(berelekamp_massey(x)[0])
    
    # Step 3 Compute mean
    a = float(M)/2.0
    b = ((((-1)**(M+1))+9.0))/36.0
    c = ((M/3.0) + (2.0/9.0))/(2**M)
    mu =  a+b-c
    
    T = list()
    for i in range(N):
        x = ((-1.0)**M) * (LC[i] - mu) + (2.0/9.0)
        T.append(x)
        
    # Step 4 Count the distribution over Ticket
    v = [0,0,0,0,0,0,0]
    for t in T:
        if t <= -2.5:
            v[0] += 1
        elif t <= -1.5:
            v[1] += 1
        elif t <= -0.5:
            v[2] += 1
        elif t <= 0.5:
            v[3] += 1
        elif t <= 1.5:
            v[4] += 1
        elif t <= 2.5:
            v[5] += 1            
        else:
            v[6] += 1

    # Step 5 Compute Chi Square Statistic
    pi = [0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]
    chi_sq = 0.0
    for i in range(K+1):
        chi_sq += ((v[i] - (N*pi[i]))**2.0)/(N*pi[i])

    # Step 6 Compute P Value
    P = ss.gammaincc((K/2.0),(chi_sq/2.0))
    # print("finished")

    success = (P >= 0.01)
    print("finished!\n")
    return [n, M, N, K, v, mu, chi_sq, P, success]
    
       
