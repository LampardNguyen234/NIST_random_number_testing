from __future__ import print_function

import math
import numpy as np
import csv
from fractions import Fraction
# from scipy.special import gamma, gammainc, gammaincc
import cmath
import random

from gamma_functions import *
from non_overlapping_template_matching_test import *
from serial_test import *

def monobit_test(input, output):

    fi = open(input, "r+")

    num_count = 0

    write_array = [0.0,0.0,0.0,0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number', 'zeroes count', 'ones count', 's', 'p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()

        for line in fi:

            num_count = num_count + 1

            n = len(line) - 1

            t = hex(int(line, 2))[2:-1].upper() #hex form of the number
            
            ones = line.count('1') #number of ones


            zeroes = line.count('0')    #number of zeros

            s = abs(ones - zeroes)  

            p = math.erfc(float(s)/(math.sqrt(float(n)) * math.sqrt(2.0))) #p-value

            if p>=0.01:
                success = "PASS"
                write_array[4] = write_array[4] + 1
            else:
                success = "FAIL"

            write_array[0] = write_array[0] + zeroes
            write_array[1] = write_array[1] + ones
            write_array[3] = write_array[3] + p

            writer.writerow({fieldnames[0]: t, fieldnames[1]: zeroes, fieldnames[2]: ones, fieldnames[3]: s, fieldnames[4]: p, fieldnames[5]:  success})

        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[1]: write_array[0]/num_count, 
            fieldnames[2]: write_array[1]/num_count, 
            fieldnames[3]: abs(write_array[0]/num_count - write_array[1]/num_count), 
            fieldnames[4]: write_array[3]/num_count, 
            fieldnames[5]: write_array[4]/num_count
            })
    return (write_array[3]/num_count, write_array[4]/num_count)

def runs_test(input, output):
    fi = open(input, "r+")

    num_count = 0

    write_array = [0.0,0.0,0.0,0.0,0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number', 'zeroes count', 'ones count', 'one_prop', 'vobs','p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()

        for line in fi:

            num_count = num_count + 1

            n = len(line) - 1

            t = hex(int(line, 2))[2:-1].upper() #hex form of the number
            
            ones = line.count('1') #number of ones

            zeroes = line.count('0')    #number of zeros

            prop = float(ones)/float(n)
    
            tau = 2.0/math.sqrt(n)

            vobs = 0.0

            if abs(prop-0.5) > tau:
                p = 0
                success = "FAIL"
            else:

                vobs = 1.0
                for i in range(n-1):
                    if line[i] != line[i+1]:
                        vobs += 1.0

            # print("  vobs ",vobs)
              
                p = math.erfc(abs(vobs - (2.0*n*prop*(1.0-prop)))/(2.0*math.sqrt(2.0*n)*prop*(1-prop) ))
                if p>=0.01:
                    success = "PASS"
                    write_array[5] = write_array[5] + 1

            write_array[0] = write_array[0] + zeroes
            write_array[1] = write_array[1] + ones
            write_array[2] = write_array[2] + prop
            write_array[3] = write_array[3] + vobs
            write_array[4] = write_array[4] + p

            writer.writerow({
                fieldnames[0]: t, 
                fieldnames[1]: zeroes, 
                fieldnames[2]: ones, 
                fieldnames[3]: prop,
                fieldnames[4]: vobs, 
                fieldnames[5]: p, 
                fieldnames[6]:  success
                })

        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[1]: write_array[0]/num_count, 
            fieldnames[2]: write_array[1]/num_count, 
            fieldnames[3]: write_array[2]/num_count, 
            fieldnames[4]: write_array[3]/num_count,
            fieldnames[5]: write_array[4]/num_count, 
            fieldnames[6]: write_array[5]/num_count
            })
    return (write_array[4]/num_count, write_array[5]/num_count)

def frequency_within_block_test(input, output):
    # Compute number of blocks M = block size. N=num of blocks
    # N = floor(n/M)
    # miniumum block size 20 bits, most blocks 100

    fi = open(input, "r+")

    num_count = 0

    write_array = [0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number','chisq','p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()

        for line in fi:

            num_count = num_count + 1

            n = len(line)-1
            
            t = hex(int(line, 2))[2:-1].upper() #hex form of the number

            M = 16
            
            N = int(math.floor(n/M))
            
            if N > 99:
                N=99
                M = int(math.floor(n/N))
        
            if n < 100:
                print("Too little data for test. Supply at least 100 bits")
                return FAIL,1.0,None
        
        
            num_of_blocks = N

            block_size = M #int(math.floor(len(bits)/num_of_blocks))
            #n = int(block_size * num_of_blocks)
        
            proportions = list()
            for i in range(num_of_blocks):
                block = line[i*(block_size):((i+1)*(block_size))]
                
                ones = block.count('1') 

                zeroes = block.count('0') 
                
                proportions.append(Fraction(ones,block_size))

            chisq = 0.0
            for prop in proportions:
                chisq += 4.0*block_size*((prop - Fraction(1,2))**2)
            
            p = gammaincc((num_of_blocks/2.0),float(chisq)/2.0)
            
            if p>=0.01:
                success = "PASS"
                write_array[1] = write_array[1] + 1
            else:
                success = "FAIL"

            write_array[0] =write_array[0] + p

            writer.writerow({
                fieldnames[0]: t, 
                fieldnames[1]: chisq, 
                fieldnames[2]: p, 
                fieldnames[3]: success
                })

        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[3]: write_array[1]/num_count,
            })

def longest_run_ones_in_a_block_test(input, output):
    fi = open(input, "r+")

    num_count = 0

    write_array = [0.0,0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number', 'chisq','p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()

        M8 =      [0.2148, 0.3672, 0.2305, 0.1875]

        M = 8
                    
        K = 3

        N = 16

        n = 256

        for line in fi:

            num_count = num_count + 1
            
            t = hex(int(line, 2))[2:-1].upper() #hex form of the number

            
                
            # Table of frequencies
            v = [0,0,0,0,0,0,0]

            for i in range(N): # over each block
                #find longest run
                block = line[i*M:((i+1)*M)] # Block i
                
                run = 0
                longest = 0
                for j in range(M): # Count the bits.
                    if block[j] == '1':
                        run += 1
                        if run > longest:
                            longest = run
                    else:
                        run = 0

                if longest <= 1:    v[0] += 1
                elif longest == 2:  v[1] += 1
                elif longest == 3:  v[2] += 1
                else:               v[3] += 1
            
            # Compute Chi-Sq
            chi_sq = 0.0
            for i in range(K+1):
                p_i = M8[i]
                upper = (v[i] - N*p_i)**2
                lower = N*p_i
                chi_sq += upper/lower

            p = gammaincc(K/2.0, chi_sq/2.0)
            
            if p>= 0.01:
                success = "PASS"
                write_array[2] = write_array[2] + 1
            else:
                success = "FAIL"

            writer.writerow({
                fieldnames[0]: t, 
                fieldnames[1]: chi_sq, 
                fieldnames[2]: p, 
                fieldnames[3]: success
                })

            write_array[1] = write_array[1] + p

        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[2]: write_array[1]/num_count,
            fieldnames[3]: write_array[2]/num_count,
            })

def dft_test(input, output, ):

    fi = open(input, "r+")

    n = 256

    num_count = 0

    T = math.sqrt(math.log(1.0/0.05)*n) # Compute upper threshold

    N0 = 0.95*n/2.0

    write_array = [0.0,0.0,0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number', 'N0', 'N1', 'd','p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()


        for line in fi:

            t = hex(int(line, 2))[2:-1].upper() #hex form of the number

            num_count = num_count + 1

            ts = list()             # Convert to +1,-1
            for i in range(n):
                if line[i] == '1':
                    ts.append(1)
                else:
                    ts.append(-1)
            ts_np = np.array(ts)

            fs = np.fft.fft(ts_np)  # Compute DFT

            mags = abs(fs)[:n/2]  #Compute magnitudes of first half of sequence


            N1 = 0.0   # Count the peaks above the upper theshold

            for mag in mags:
                if mag < T:
                    N1 += 1.0
            d = (N1 - N0)/math.sqrt((n*0.95*0.05)/4) # Compute the P value

            p = math.erfc(abs(d)/math.sqrt(2))

            if p>=0.01:
                success = "PASS"
                write_array[3] = write_array[3] + 1
            else:
                success = "FAIL"

            writer.writerow({
                fieldnames[0]: t, 
                fieldnames[1]: N0, 
                fieldnames[2]: N1, 
                fieldnames[3]: d,
                fieldnames[4]: p, 
                fieldnames[5]: success
                })
        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[5]: write_array[3]/num_count,
            })

def non_overlapping_test(input, output):
    fi = open(input, "r+")

    n = 256

    num_count = 0

    write_array = [0.0,0.0,0.0,0.0]

    with open(output, mode="w") as fo:
        fieldnames = ['number', 'mu', 'sigma', 'chi_sq','p-value', 'success']

        writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()


        for line in fi:

            t = hex(int(line, 2))[2:-1].upper() #hex form of the number

            num_count = num_count + 1

            mu, sigma, chisq, p, success = non_overlapping_template_matching_test(line, n)

            writer.writerow({
                fieldnames[0]: t, 
                fieldnames[1]: mu, 
                fieldnames[2]: sigma, 
                fieldnames[3]: chisq,
                fieldnames[4]: p, 
                fieldnames[5]: success
                })

            if p>=0.01:
                write_array[3] += 1

        writer.writerow({
            fieldnames[0]: "Average",
            fieldnames[5]: write_array[3]/num_count,
            })


def serial_test1(input, output):
	
	fi = open(input, "r+")

	num_count = 0
	
	write_array = [0.0,0.0,0.0,0.0]

   	with open(output, mode="w") as fo:

   		fieldnames = ['number', 'mu', 'sigma', 'chi_sq','p-value', 'success']

		writer = csv.DictWriter(fo, fieldnames=fieldnames)

		writer.writeheader()

		for line in fi:

			t = hex(int(line, 2))[2:-1].upper() #hex form of the number

    	   	num_count += 1

    	   	result = serial_test(line)



def main():
    input = "input.txt"

    output_monobit = "result_monobit_test.csv"

    output_runs = "result_runs_test.csv"

    output_fwb = "result_fwb_test.csv"

    output_longest_run_in_a_block = "result_longest_run_in_a_block_test.csv"

    output_dft = "result_dft_test.csv"

    output_nonoverlapping = "result_nonoverlapping_test.csv"

    output_serial = "result_serial_test.csv"


    # monobit_test(input,output_monobit)

    # runs_test(input,output_runs)

    # frequency_within_block_test(input, output_fwb)

    # longest_run_ones_in_a_block_test(input, output_longest_run_in_a_block)

    # dft_test(input, output_dft)

    # non_overlapping_test(input, output_nonoverlapping)

    serial_test1(input,output_serial)


if __name__ == "__main__":
    main()