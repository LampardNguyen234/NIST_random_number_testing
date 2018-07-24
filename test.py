import csv

testlist = [
        'monobit_test',
        'frequency_within_block_test',
        'runs_test',
        'longest_run_ones_in_a_block_test',
        'binary_matrix_rank_test',
        'dft_test',
        'non_overlapping_template_matching_test',
        'overlapping_template_matching_test',
        'maurers_universal_test',
        'linear_complexity_test',
        'serial_test',
        'approximate_entropy_test',
        'cumulative_sums_test',
        'random_excursion_test',
        'random_excursion_variant_test'
        ]

def main():
    input = "input.txt"

    NUM_TEST = 15

    fieldnames = [None] * NUM_TEST    #fieldnames according to csv file

    fieldnames[0] = ['number', 'zeroes count', 'ones count', 's', 'p-value', 'success']
    fieldnames[1] = ['number','chisq','p-value', 'success']
    fieldnames[2] = ['number', 'zeroes count', 'ones count', 'one_prop', 'vobs','p-value', 'success']
    fieldnames[3] = ['number', 'chisq','p-value', 'success']
    fieldnames[5] = ['number', 'N0', 'N1', 'd','p-value', 'success']
    fieldnames[6] = ['number', 'mu', 'sigma', 'chi_sq','p-value', 'success']
    fieldnames[10] = ['number', 'psi_sq_m', 'psi_sq_mm1', 'psi_sq_mm2', 'delta1', 'delta2', 'p1', 'p2', 'success']
    fieldnames[11] = ['number', 'appen_m', 'chi_sqp', 'p', 'success']
    fieldnames[12] = ['number', 'p_forward', 'p_backward', 'success']
    fieldnames[13] = ['number', 'J', 'success']
    
    fo = [None]*NUM_TEST #file out

    output = [None]*NUM_TEST  #output file name
    writer = [None]*NUM_TEST  #writers of csv file
 
    source = [None]*NUM_TEST

    for i in range(NUM_TEST):
        if i<9:
            output[i] = "results/result_0" + str(i+1) + "_"+ testlist[i] + ".csv"
        else:
            output[i] = "results/result_" + str(i+1) + "_" + testlist[i] + ".csv"

        if fieldnames[i] != None:
            fo[i] = open(output[i], mode="w")
            writer[i] = csv.DictWriter(fo[i], fieldnames=fieldnames[i])
            writer[i].writeheader()

    fi = open(input, "r+") # input file

    result = [None]*NUM_TEST


    for line in fi:
        for i in range(NUM_TEST):
            if fieldnames[i] != None:
                # Get test file
                if i<9:
                    m = __import__("0"+ str(i+1) + "_" + testlist[i])
                else:
                    m = __import__(str(i+1) + "_" + testlist[i])

                x = m.test(line[:-1], 256)

                print(testlist[i] + ": result = " + str(x[len(x)-1]))
        exit()
        # for i in range(NUM_TEST):




if __name__ == "__main__":
    main()