#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

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
    outputfile = "../../solutions/%s/%s-thop/%s_%02d_%s_%s_%02d_%02d.thop.sol" % ("ils", tsp_base, tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition+1)
    os.system("./ilsthop --inputfile %s --seed %d --time %.1f --outputfile %s" % (inputfile, 
                                                                                  random_seeds[repetition], 
                                                                                  math.ceil((int(''.join(filter(lambda x: x.isdigit(), tsp_base))) - 2) * number_of_items_per_city / 10.0), 
                                                                                  outputfile
                                                                                 ))

if __name__ == "__main__":

    tsp_base = ["eil51", "pr107", "a280", "dsj1000", ]
    number_of_items_per_city = [1, 3, 5, 10, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, 5, 10, ]
    maximum_travel_time = [1, 2, 3, ]
    number_of_runs = 30

    os.system("make clean")
    os.system("make ils")

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
