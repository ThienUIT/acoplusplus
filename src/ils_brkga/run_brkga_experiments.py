#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

parameter_configurations = {
"eil51_01_bsc": {"--p": 4687, "--pe": 0.36, "--pm": 0.1, "--rhoe": 0.66, "--lsfreq": 938},
"eil51_01_unc": {"--p": 3427, "--pe": 0.55, "--pm": 0.11, "--rhoe": 0.72, "--lsfreq": 594},
"eil51_01_usw": {"--p": 3930, "--pe": 0.31, "--pm": 0.16, "--rhoe": 0.58, "--lsfreq": 756},
"eil51_03_bsc": {"--p": 756, "--pe": 0.24, "--pm": 0.05, "--rhoe": 0.54, "--lsfreq": 353},
"eil51_03_unc": {"--p": 750, "--pe": 0.28, "--pm": 0.07, "--rhoe": 0.62, "--lsfreq": 264},
"eil51_03_usw": {"--p": 791, "--pe": 0.32, "--pm": 0.07, "--rhoe": 0.61, "--lsfreq": 415},
"eil51_05_bsc": {"--p": 510, "--pe": 0.26, "--pm": 0.03, "--rhoe": 0.55, "--lsfreq": 593},
"eil51_05_unc": {"--p": 376, "--pe": 0.28, "--pm": 0.08, "--rhoe": 0.56, "--lsfreq": 849},
"eil51_05_usw": {"--p": 291, "--pe": 0.22, "--pm": 0.05, "--rhoe": 0.53, "--lsfreq": 709},
"eil51_10_bsc": {"--p": 491, "--pe": 0.13, "--pm": 0.06, "--rhoe": 0.69, "--lsfreq": 250},
"eil51_10_unc": {"--p": 171, "--pe": 0.27, "--pm": 0.03, "--rhoe": 0.63, "--lsfreq": 701},
"eil51_10_usw": {"--p": 210, "--pe": 0.28, "--pm": 0.03, "--rhoe": 0.7, "--lsfreq": 417},
"pr107_01_bsc": {"--p": 4986, "--pe": 0.27, "--pm": 0.09, "--rhoe": 0.61, "--lsfreq": 384},
"pr107_01_unc": {"--p": 4932, "--pe": 0.34, "--pm": 0.05, "--rhoe": 0.53, "--lsfreq": 224},
"pr107_01_usw": {"--p": 4855, "--pe": 0.43, "--pm": 0.06, "--rhoe": 0.64, "--lsfreq": 297},
"pr107_03_bsc": {"--p": 619, "--pe": 0.3, "--pm": 0.07, "--rhoe": 0.56, "--lsfreq": 677},
"pr107_03_unc": {"--p": 707, "--pe": 0.17, "--pm": 0.2, "--rhoe": 0.63, "--lsfreq": 499},
"pr107_03_usw": {"--p": 407, "--pe": 0.31, "--pm": 0.09, "--rhoe": 0.58, "--lsfreq": 684},
"pr107_05_bsc": {"--p": 260, "--pe": 0.29, "--pm": 0.08, "--rhoe": 0.59, "--lsfreq": 485},
"pr107_05_unc": {"--p": 409, "--pe": 0.22, "--pm": 0.09, "--rhoe": 0.59, "--lsfreq": 469},
"pr107_05_usw": {"--p": 271, "--pe": 0.44, "--pm": 0.04, "--rhoe": 0.64, "--lsfreq": 827},
"pr107_10_bsc": {"--p": 173, "--pe": 0.24, "--pm": 0.18, "--rhoe": 0.66, "--lsfreq": 632},
"pr107_10_unc": {"--p": 273, "--pe": 0.31, "--pm": 0.04, "--rhoe": 0.62, "--lsfreq": 549},
"pr107_10_usw": {"--p": 205, "--pe": 0.26, "--pm": 0.11, "--rhoe": 0.67, "--lsfreq": 614},
"a280_01_bsc": {"--p": 3407, "--pe": 0.43, "--pm": 0.03, "--rhoe": 0.63, "--lsfreq": 396},
"a280_01_unc": {"--p": 4631, "--pe": 0.36, "--pm": 0.01, "--rhoe": 0.59, "--lsfreq": 219},
"a280_01_usw": {"--p": 4193, "--pe": 0.31, "--pm": 0.03, "--rhoe": 0.62, "--lsfreq": 211},
"a280_03_bsc": {"--p": 558, "--pe": 0.45, "--pm": 0.01, "--rhoe": 0.63, "--lsfreq": 623},
"a280_03_unc": {"--p": 313, "--pe": 0.54, "--pm": 0.03, "--rhoe": 0.67, "--lsfreq": 862},
"a280_03_usw": {"--p": 382, "--pe": 0.49, "--pm": 0.01, "--rhoe": 0.64, "--lsfreq": 985},
"a280_05_bsc": {"--p": 222, "--pe": 0.39, "--pm": 0.03, "--rhoe": 0.56, "--lsfreq": 794},
"a280_05_unc": {"--p": 318, "--pe": 0.26, "--pm": 0.03, "--rhoe": 0.61, "--lsfreq": 704},
"a280_05_usw": {"--p": 274, "--pe": 0.33, "--pm": 0.04, "--rhoe": 0.67, "--lsfreq": 764},
"a280_10_bsc": {"--p": 165, "--pe": 0.27, "--pm": 0.01, "--rhoe": 0.67, "--lsfreq": 693},
"a280_10_unc": {"--p": 188, "--pe": 0.21, "--pm": 0.03, "--rhoe": 0.61, "--lsfreq": 821},
"a280_10_usw": {"--p": 210, "--pe": 0.37, "--pm": 0.01, "--rhoe": 0.65, "--lsfreq": 789},
"dsj1000_01_bsc": {"--p": 4233, "--pe": 0.31, "--pm": 0, "--rhoe": 0.56, "--lsfreq": 569},
"dsj1000_01_unc": {"--p": 3663, "--pe": 0.31, "--pm": 0.01, "--rhoe": 0.55, "--lsfreq": 646},
"dsj1000_01_usw": {"--p": 3300, "--pe": 0.44, "--pm": 0.01, "--rhoe": 0.57, "--lsfreq": 604},
"dsj1000_03_bsc": {"--p": 182, "--pe": 0.43, "--pm": 0, "--rhoe": 0.66, "--lsfreq": 587},
"dsj1000_03_unc": {"--p": 241, "--pe": 0.16, "--pm": 0.03, "--rhoe": 0.55, "--lsfreq": 340},
"dsj1000_03_usw": {"--p": 136, "--pe": 0.4, "--pm": 0, "--rhoe": 0.58, "--lsfreq": 565},
"dsj1000_05_bsc": {"--p": 142, "--pe": 0.24, "--pm": 0.02, "--rhoe": 0.61, "--lsfreq": 658},
"dsj1000_05_unc": {"--p": 117, "--pe": 0.41, "--pm": 0, "--rhoe": 0.55, "--lsfreq": 737},
"dsj1000_05_usw": {"--p": 163, "--pe": 0.32, "--pm": 0.01, "--rhoe": 0.64, "--lsfreq": 579},
"dsj1000_10_bsc": {"--p": 86, "--pe": 0.26, "--pm": 0.02, "--rhoe": 0.69, "--lsfreq": 316},
"dsj1000_10_unc": {"--p": 112, "--pe": 0.19, "--pm": 0.01, "--rhoe": 0.65, "--lsfreq": 271},
"dsj1000_10_usw": {"--p": 71, "--pe": 0.36, "--pm": 0.01, "--rhoe": 0.62, "--lsfreq": 715},
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

def launcher(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition):
    if knapsack_size != "inf": knapsack_size = "%02d" % (knapsack_size, )
    inputfile = "../../instances/%s-thop/%s_%02d_%s_%s_%02d.thop" % (tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time)
    outputfile = "../../solutions/%s/%s-thop/%s_%02d_%s_%s_%02d_%02d.thop.sol" % ("brkga", tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition+1)
    parameter_configuration_key = "%s_%02d_%s" % (tsp_base, number_of_items_per_city, knapsack_type)

    os.system("./brkgathop --inputfile %s --seed %d --time %.1f --outputfile %s %s" % (inputfile, 
                                                                                       random_seeds[repetition], 
                                                                                       math.ceil((int(''.join(filter(lambda x: x.isdigit(), tsp_base))) - 2) * number_of_items_per_city / 10.0), 
                                                                                       outputfile,
                                                                                       ' '.join("%s %s" % (k, v) for k, v in parameter_configurations[parameter_configuration_key].items())
                                                                                      ))

if __name__ == "__main__":

    tsp_base = ["eil51", "pr107", "a280", "dsj1000", ]
    number_of_items_per_city = [1, 3, 5, 10, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, 5, 10, ]
    maximum_travel_time = [1, 2, 3, ]
    number_of_runs = 30

    os.system("make clean")
    os.system("make brkga")

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
