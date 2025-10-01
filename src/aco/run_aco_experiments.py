#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

parameter_configurations = {
"eil51_01_bsc": {"--ants": 100, "--alpha": 0.96, "--beta": 2.72, "--rho": 0.69, "--ptries": 3},
"eil51_01_unc": {"--ants": 100, "--alpha": 0.97, "--beta": 3.82, "--rho": 0.87, "--ptries": 2},
"eil51_01_usw": {"--ants": 200, "--alpha": 0.79, "--beta": 4.46, "--rho": 0.92, "--ptries": 1},
"eil51_03_bsc": {"--ants": 100, "--alpha": 0.71, "--beta": 4.92, "--rho": 0.75, "--ptries": 3},
"eil51_03_unc": {"--ants": 100, "--alpha": 0.8, "--beta": 3.89, "--rho": 0.78, "--ptries": 1},
"eil51_03_usw": {"--ants": 100, "--alpha": 0.86, "--beta": 3.64, "--rho": 0.82, "--ptries": 1},
"eil51_05_bsc": {"--ants": 100, "--alpha": 0.83, "--beta": 4.23, "--rho": 0.86, "--ptries": 2},
"eil51_05_unc": {"--ants": 100, "--alpha": 0.9, "--beta": 4.83, "--rho": 0.88, "--ptries": 2},
"eil51_05_usw": {"--ants": 200, "--alpha": 0.89, "--beta": 2.69, "--rho": 0.78, "--ptries": 1},
"eil51_10_bsc": {"--ants": 100, "--alpha": 0.95, "--beta": 4.81, "--rho": 0.91, "--ptries": 2},
"eil51_10_unc": {"--ants": 100, "--alpha": 0.94, "--beta": 4.66, "--rho": 0.91, "--ptries": 2},
"eil51_10_usw": {"--ants": 50, "--alpha": 0.88, "--beta": 5.43, "--rho": 0.74, "--ptries": 1},
"pr107_01_bsc": {"--ants": 1000, "--alpha": 0.86, "--beta": 2.37, "--rho": 0.65, "--ptries": 2},
"pr107_01_unc": {"--ants": 100, "--alpha": 0.8, "--beta": 3.18, "--rho": 0.59, "--ptries": 4},
"pr107_01_usw": {"--ants": 100, "--alpha": 0.97, "--beta": 3.56, "--rho": 0.9, "--ptries": 1},
"pr107_03_bsc": {"--ants": 200, "--alpha": 0.8, "--beta": 3.48, "--rho": 0.32, "--ptries": 3},
"pr107_03_unc": {"--ants": 200, "--alpha": 0.81, "--beta": 2.97, "--rho": 0.58, "--ptries": 2},
"pr107_03_usw": {"--ants": 200, "--alpha": 0.86, "--beta": 3.27, "--rho": 0.49, "--ptries": 2},
"pr107_05_bsc": {"--ants": 200, "--alpha": 0.81, "--beta": 2.41, "--rho": 0.37, "--ptries": 2},
"pr107_05_unc": {"--ants": 200, "--alpha": 0.91, "--beta": 4.72, "--rho": 0.73, "--ptries": 3},
"pr107_05_usw": {"--ants": 500, "--alpha": 0.88, "--beta": 3.12, "--rho": 0.9, "--ptries": 1},
"pr107_10_bsc": {"--ants": 100, "--alpha": 0.91, "--beta": 4.03, "--rho": 0.76, "--ptries": 3},
"pr107_10_unc": {"--ants": 500, "--alpha": 0.82, "--beta": 2.82, "--rho": 0.64, "--ptries": 2},
"pr107_10_usw": {"--ants": 200, "--alpha": 0.85, "--beta": 3.1, "--rho": 0.56, "--ptries": 1},
"a280_01_bsc": {"--ants": 200, "--alpha": 0.71, "--beta": 4.34, "--rho": 0.7, "--ptries": 2},
"a280_01_unc": {"--ants": 100, "--alpha": 0.79, "--beta": 5.28, "--rho": 0.75, "--ptries": 1},
"a280_01_usw": {"--ants": 200, "--alpha": 0.74, "--beta": 5.27, "--rho": 0.54, "--ptries": 1},
"a280_03_bsc": {"--ants": 200, "--alpha": 0.75, "--beta": 7.06, "--rho": 0.51, "--ptries": 2},
"a280_03_unc": {"--ants": 500, "--alpha": 0.95, "--beta": 5.26, "--rho": 0.51, "--ptries": 2},
"a280_03_usw": {"--ants": 200, "--alpha": 0.66, "--beta": 6.29, "--rho": 0.58, "--ptries": 1},
"a280_05_bsc": {"--ants": 100, "--alpha": 0.72, "--beta": 9.15, "--rho": 0.9, "--ptries": 2},
"a280_05_unc": {"--ants": 100, "--alpha": 0.72, "--beta": 5.72, "--rho": 0.71, "--ptries": 3},
"a280_05_usw": {"--ants": 50, "--alpha": 0.74, "--beta": 6.46, "--rho": 0.73, "--ptries": 1},
"a280_10_bsc": {"--ants": 100, "--alpha": 0.78, "--beta": 5.02, "--rho": 0.6, "--ptries": 1},
"a280_10_unc": {"--ants": 20, "--alpha": 0.99, "--beta": 4.93, "--rho": 0.4, "--ptries": 1},
"a280_10_usw": {"--ants": 100, "--alpha": 0.84, "--beta": 7.52, "--rho": 0.64, "--ptries": 1},
"dsj1000_01_bsc": {"--ants": 200, "--alpha": 0.89, "--beta": 7.27, "--rho": 0.63, "--ptries": 1},
"dsj1000_01_unc": {"--ants": 100, "--alpha": 0.89, "--beta": 7.76, "--rho": 0.54, "--ptries": 2},
"dsj1000_01_usw": {"--ants": 50, "--alpha": 0.9, "--beta": 5.06, "--rho": 0.53, "--ptries": 3},
"dsj1000_03_bsc": {"--ants": 100, "--alpha": 0.94, "--beta": 4.79, "--rho": 0.34, "--ptries": 1},
"dsj1000_03_unc": {"--ants": 20, "--alpha": 0.92, "--beta": 6.92, "--rho": 0.39, "--ptries": 5},
"dsj1000_03_usw": {"--ants": 100, "--alpha": 0.89, "--beta": 6.45, "--rho": 0.69, "--ptries": 1},
"dsj1000_05_bsc": {"--ants": 50, "--alpha": 0.9, "--beta": 7.15, "--rho": 0.45, "--ptries": 1},
"dsj1000_05_unc": {"--ants": 200, "--alpha": 1.69, "--beta": 7.44, "--rho": 0.27, "--ptries": 4},
"dsj1000_05_usw": {"--ants": 100, "--alpha": 0.87, "--beta": 6.15, "--rho": 0.47, "--ptries": 1},
"dsj1000_10_bsc": {"--ants": 50, "--alpha": 1.13, "--beta": 8.3, "--rho": 0.24, "--ptries": 1},
"dsj1000_10_unc": {"--ants": 200, "--alpha": 1.02, "--beta": 7.85, "--rho": 0.39, "--ptries": 2},
"dsj1000_10_usw": {"--ants": 200, "--alpha": 0.9, "--beta": 5.83, "--rho": 0.54, "--ptries": 1},
}

random_seeds = [ 269070,  99470, 126489, 644764, 547617, 642580,  73456, 462018, 858990, 756112, 
                 701531, 342080, 613485, 131654, 886148, 909040, 146518, 782904,   3075, 974703, 
                 170425, 531298, 253045, 488197, 394197, 519912, 606939, 480271, 117561, 900952, 
                 968235, 345118, 750253, 420440, 761205, 130467, 928803, 768798, 640300, 871462, 
                 639622,  90614, 187822, 594363, 193911, 846042, 680779, 344008, 759862, 661168, 
                 223420, 959508,  62985, 349296, 910428, 964420, 422964, 384194, 985214,  57575, 
                 639619,  90505, 435236, 465842, 102567, 189997, 741017, 611828, 699223, 335142, 
                  52119,  49256, 324523, 348215, 651525, 517999, 830566, 958538, 880422, 390645, 
                 148265, 807740, 934464, 524847, 408760, 668587, 257030, 751580,  90477, 594476, 
                 571216, 306614, 308010, 661191, 890429, 425031,  69108, 435783,  17725, 335928, ]

def launcher(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition, runtime_factor="1x"):
    if knapsack_size != "inf": knapsack_size = "%02d" % (knapsack_size, )
    inputfile = "../../instances/%s-thop/%s_%02d_%s_%s_%02d.thop" % (tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time)
    outputfile = "../../solutions/aco/%s-thop/%s_%02d_%s_%s_%02d_%02d.thop.sol" % (tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition+1)
    parameter_configuration_key = "%s_%02d_%s" % (tsp_base, number_of_items_per_city, knapsack_type) 
    os.system("./acothop --mmas --tries 1 --seed %d --time %.1f --inputfile %s --outputfile %s %s --log" % (random_seeds[repetition], \
                                                                                                        float(runtime_factor.replace('x','')) * math.ceil((int(''.join(filter(lambda x: x.isdigit(), tsp_base))) - 2) * number_of_items_per_city / 10.0), \
                                                                                                        inputfile, outputfile, ' '.join("%s %s" % (k, v) for k, v in parameter_configurations[parameter_configuration_key].items())))

if __name__ == "__main__":

    tsp_base = ["eil51", "pr107", "a280", "dsj1000", ]
    number_of_items_per_city = [1, 3, 5, 10, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, 5, 10, ]
    maximum_travel_time = [1, 2, 3, ]
    number_of_runs = 30
    
    os.system("make clean")
    os.system("make")

    pool = multiprocessing.Pool(processes=max(1, multiprocessing.cpu_count() - 2))
    
    for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
        _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
        for repetition in range(number_of_runs):
            pool.apply_async(launcher, args=(_tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time, repetition))         

    number_of_items_per_city = [1, ]
    knapsack_size = ["inf", ]

    for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
        _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
        for repetition in range(number_of_runs):
            pool.apply_async(launcher, args=(_tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time, repetition)) 

    pool.close()
    pool.join()
