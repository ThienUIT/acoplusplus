#!/usr/bin/python
# -*- coding: utf-8 -*-

import itertools
import os
import multiprocessing
import argparse
import math

EPS = 10e-4

def create_op_instance(thop_input_file, op_output_file):

    file = open(thop_input_file, 'r')
    lines = file.readlines()
    file.close()
  
    for i in range(len(lines)):
        if "CAPACITY OF KNAPSACK" in lines[i]:
            lines[i] = "CAPACITY OF KNAPSACK: 987654321"
        if "MIN SPEED" in lines[i]:
            lines[i] = "MIN SPEED: 1.00"
        if "MAX SPEED" in lines[i]:
            lines[i] = "MAX SPEED: 1.00"
            
    file = open(op_output_file, 'w')
    for line in lines:
        file.write(line.strip() + "\n")
    file.close()

if __name__ == "__main__":

    tsp_base = ["eil51", "pr107", "a280", "dsj1000", ]
    number_of_items_per_city = [1, ]
    knapsack_type = ["bsc", "unc", "usw", ]
    knapsack_size = [1, ]
    maximum_travel_time = [1, 2, 3, ]

    for _product in itertools.product(tsp_base, number_of_items_per_city, knapsack_type, knapsack_size, maximum_travel_time):
        _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time = _product
        input_file = "%s-thop/%s_%02d_%s_%02d_%02d.thop" % (_tsp_base, _tsp_base, _number_of_items_per_city, _knapsack_type, _knapsack_size, _maximum_travel_time)
        output_file = "%s-thop/%s_%02d_%s_inf_%02d.thop" % (_tsp_base, _tsp_base, _number_of_items_per_city, _knapsack_type, _maximum_travel_time)        
        create_op_instance(input_file, output_file)        
