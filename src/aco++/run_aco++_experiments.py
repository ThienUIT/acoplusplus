#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

parameter_configurations = {
"eil51_01_bsc": {"--ants": " 500", "--alpha": "0.66", "--beta": "4.42", "--rho": "0.90", "--ptries": "2", "--localsearch": "0"}, 
"eil51_01_unc": {"--ants": "1000", "--alpha": "0.69", "--beta": "4.29", "--rho": "0.75", "--ptries": "2", "--localsearch": "1"}, 
"eil51_01_usw": {"--ants": " 100", "--alpha": "0.92", "--beta": "2.90", "--rho": "0.19", "--ptries": "1", "--localsearch": "1"}, 
"eil51_03_bsc": {"--ants": " 200", "--alpha": "0.78", "--beta": "3.65", "--rho": "0.92", "--ptries": "1", "--localsearch": "1"}, 
"eil51_03_unc": {"--ants": " 200", "--alpha": "0.89", "--beta": "3.01", "--rho": "0.14", "--ptries": "1", "--localsearch": "0"}, 
"eil51_03_usw": {"--ants": " 100", "--alpha": "0.80", "--beta": "5.26", "--rho": "0.84", "--ptries": "1", "--localsearch": "1"}, 
"eil51_05_bsc": {"--ants": " 500", "--alpha": "0.79", "--beta": "3.12", "--rho": "0.35", "--ptries": "5", "--localsearch": "1"}, 
"eil51_05_unc": {"--ants": " 100", "--alpha": "0.88", "--beta": "3.65", "--rho": "0.78", "--ptries": "1", "--localsearch": "1"}, 
"eil51_05_usw": {"--ants": " 500", "--alpha": "0.99", "--beta": "1.38", "--rho": "0.34", "--ptries": "2", "--localsearch": "1"}, 
"eil51_10_bsc": {"--ants": " 200", "--alpha": "0.88", "--beta": "5.09", "--rho": "0.86", "--ptries": "2", "--localsearch": "1"}, 
"eil51_10_unc": {"--ants": " 100", "--alpha": "0.87", "--beta": "4.32", "--rho": "0.94", "--ptries": "1", "--localsearch": "1"}, 
"eil51_10_usw": {"--ants": " 100", "--alpha": "0.86", "--beta": "3.89", "--rho": "0.59", "--ptries": "1", "--localsearch": "1"}, 
"pr107_01_bsc": {"--ants": " 200", "--alpha": "1.64", "--beta": "1.67", "--rho": "0.18", "--ptries": "4", "--localsearch": "1"}, 
"pr107_01_unc": {"--ants": "1000", "--alpha": "2.87", "--beta": "3.92", "--rho": "0.34", "--ptries": "5", "--localsearch": "1"}, 
"pr107_01_usw": {"--ants": " 500", "--alpha": "0.93", "--beta": "2.59", "--rho": "0.60", "--ptries": "3", "--localsearch": "1"}, 
"pr107_03_bsc": {"--ants": " 200", "--alpha": "1.08", "--beta": "2.90", "--rho": "0.36", "--ptries": "4", "--localsearch": "1"}, 
"pr107_03_unc": {"--ants": " 500", "--alpha": "1.08", "--beta": "5.43", "--rho": "0.46", "--ptries": "3", "--localsearch": "1"}, 
"pr107_03_usw": {"--ants": " 500", "--alpha": "0.92", "--beta": "2.35", "--rho": "0.68", "--ptries": "3", "--localsearch": "1"}, 
"pr107_05_bsc": {"--ants": "1000", "--alpha": "1.38", "--beta": "2.60", "--rho": "0.29", "--ptries": "2", "--localsearch": "1"}, 
"pr107_05_unc": {"--ants": "1000", "--alpha": "3.54", "--beta": "4.54", "--rho": "0.20", "--ptries": "3", "--localsearch": "1"}, 
"pr107_05_usw": {"--ants": " 500", "--alpha": "1.00", "--beta": "1.73", "--rho": "0.51", "--ptries": "2", "--localsearch": "1"}, 
"pr107_10_bsc": {"--ants": " 100", "--alpha": "0.87", "--beta": "2.97", "--rho": "0.74", "--ptries": "5", "--localsearch": "1"}, 
"pr107_10_unc": {"--ants": " 200", "--alpha": "0.91", "--beta": "2.80", "--rho": "0.63", "--ptries": "4", "--localsearch": "1"}, 
"pr107_10_usw": {"--ants": "  20", "--alpha": "0.90", "--beta": "3.20", "--rho": "0.12", "--ptries": "3", "--localsearch": "1"}, 
"a280_01_bsc": {"--ants": " 200", "--alpha": "0.71", "--beta": "6.93", "--rho": "0.42", "--ptries": "2", "--localsearch": "0"}, 
"a280_01_unc": {"--ants": " 200", "--alpha": "1.16", "--beta": "1.79", "--rho": "0.42", "--ptries": "3", "--localsearch": "1"}, 
"a280_01_usw": {"--ants": " 200", "--alpha": "0.84", "--beta": "3.83", "--rho": "0.39", "--ptries": "2", "--localsearch": "0"}, 
"a280_03_bsc": {"--ants": " 100", "--alpha": "0.67", "--beta": "7.12", "--rho": "0.53", "--ptries": "1", "--localsearch": "0"}, 
"a280_03_unc": {"--ants": " 200", "--alpha": "0.75", "--beta": "4.10", "--rho": "0.74", "--ptries": "5", "--localsearch": "1"}, 
"a280_03_usw": {"--ants": "  20", "--alpha": "0.80", "--beta": "5.83", "--rho": "0.25", "--ptries": "2", "--localsearch": "1"}, 
"a280_05_bsc": {"--ants": "1000", "--alpha": "0.66", "--beta": "6.46", "--rho": "0.31", "--ptries": "1", "--localsearch": "0"}, 
"a280_05_unc": {"--ants": " 200", "--alpha": "0.84", "--beta": "5.89", "--rho": "0.24", "--ptries": "2", "--localsearch": "1"}, 
"a280_05_usw": {"--ants": " 200", "--alpha": "0.86", "--beta": "7.94", "--rho": "0.42", "--ptries": "2", "--localsearch": "1"}, 
"a280_10_bsc": {"--ants": " 200", "--alpha": "0.77", "--beta": "5.71", "--rho": "0.27", "--ptries": "2", "--localsearch": "1"}, 
"a280_10_unc": {"--ants": " 200", "--alpha": "0.86", "--beta": "6.12", "--rho": "0.30", "--ptries": "3", "--localsearch": "1"}, 
"a280_10_usw": {"--ants": " 200", "--alpha": "0.77", "--beta": "6.30", "--rho": "0.35", "--ptries": "2", "--localsearch": "1"}, 
"dsj1000_01_bsc": {"--ants": " 200", "--alpha": "4.30", "--beta": "9.18", "--rho": "0.27", "--ptries": "4", "--localsearch": "1"}, 
"dsj1000_01_unc": {"--ants": "  50", "--alpha": "8.33", "--beta": "6.39", "--rho": "0.10", "--ptries": "5", "--localsearch": "1"}, 
"dsj1000_01_usw": {"--ants": " 100", "--alpha": "2.49", "--beta": "6.19", "--rho": "0.37", "--ptries": "5", "--localsearch": "1"}, 
"dsj1000_03_bsc": {"--ants": "  50", "--alpha": "0.99", "--beta": "5.16", "--rho": "0.81", "--ptries": "4", "--localsearch": "1"}, 
"dsj1000_03_unc": {"--ants": " 100", "--alpha": "2.71", "--beta": "8.81", "--rho": "0.25", "--ptries": "2", "--localsearch": "1"}, 
"dsj1000_03_usw": {"--ants": " 100", "--alpha": "0.90", "--beta": "5.03", "--rho": "0.84", "--ptries": "5", "--localsearch": "1"}, 
"dsj1000_05_bsc": {"--ants": "  20", "--alpha": "0.88", "--beta": "6.10", "--rho": "0.69", "--ptries": "3", "--localsearch": "1"}, 
"dsj1000_05_unc": {"--ants": " 100", "--alpha": "5.93", "--beta": "8.50", "--rho": "0.12", "--ptries": "5", "--localsearch": "1"}, 
"dsj1000_05_usw": {"--ants": " 100", "--alpha": "0.88", "--beta": "6.18", "--rho": "0.83", "--ptries": "2", "--localsearch": "1"}, 
"dsj1000_10_bsc": {"--ants": " 500", "--alpha": "3.20", "--beta": "7.89", "--rho": "0.21", "--ptries": "1", "--localsearch": "1"}, 
"dsj1000_10_unc": {"--ants": " 500", "--alpha": "2.94", "--beta": "7.96", "--rho": "0.26", "--ptries": "1", "--localsearch": "1"}, 
"dsj1000_10_usw": {"--ants": " 500", "--alpha": "4.04", "--beta": "7.82", "--rho": "0.39", "--ptries": "2", "--localsearch": "1"},
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

def launcher(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time, repetition, runtime_factor="1t"):
    if knapsack_size != "inf":
        knapsack_size = "%02d" % (knapsack_size, )

    # Build input and output file paths
    inputfile = os.path.normpath(os.path.join(
        os.path.dirname(__file__), "..", "..", "instances", f"{tsp_base}-thop",
        f"{tsp_base}_{number_of_items_per_city:02d}_{knapsack_type}_{knapsack_size}_{maximum_travel_time:02d}.thop"
    ))
    output_dir = os.path.normpath(os.path.join(
        os.path.dirname(__file__), "..", "..", "solutions", "aco++", f"{tsp_base}-thop"
    ))
    outputfile = os.path.join(
        output_dir,
        f"{tsp_base}_{number_of_items_per_city:02d}_{knapsack_type}_{knapsack_size}_{maximum_travel_time:02d}_{repetition+1:02d}.thop.sol"
    )

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    parameter_configuration_key = f"{tsp_base}_{number_of_items_per_city:02d}_{knapsack_type}"

    # Use the correct executable name for Windows
    exe_name = "acothop.exe" if os.name == "nt" else "./acothop"

    cmd = (
        f"{exe_name} --mmas --tries 1 --seed {random_seeds[repetition]} "
        f"--time {float(runtime_factor.replace('t','')) * math.ceil((int(''.join(filter(str.isdigit, tsp_base))) - 2) * number_of_items_per_city / 10.0):.1f} "
        f'--inputfile "{inputfile}" --outputfile "{outputfile}" '
        f'{" ".join(f"{k} {v}" for k, v in parameter_configurations[parameter_configuration_key].items())} --log'
    )
    os.system(cmd)
    
if __name__ == "__main__":

    tsp_base = ["eil51" ]
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
