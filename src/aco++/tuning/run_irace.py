import numpy as np
import argparse
import os
import math

scenario =  "###################################################### -*- mode: r -*- #####\n"\
            "## Scenario setup for Iterated Race (iRace).\n"\
            "############################################################################\n"\
            "## To use the default value of a parameter of iRace, simply do not set\n"\
            "## the parameter (comment it out in this file, and do not give any\n"\
            "## value on the command line).\n"\
            "\n"\
            "## File that contains the description of the parameters.\n"\
            "parameterFile = \"./parameters_aco++.txt\"\n"\
            "\n"\
            "## Directory where the programs will be run.\n"\
            "execDir = \"./\"\n"\
            "\n"\
            "## Directory where tuning instances are located, either absolute path or\n"\
            "## relative to current directory.\n"\
            "trainInstancesDir = \"./train_instances\"\n"\
            "\n"\
            "## File that contains a list of logical expressions that cannot be TRUE\n"\
            "## for any evaluated configuration. If empty or NULL, do not use forbidden\n"\
            "## expressions.\n"\
            "forbiddenFile = \"./forbidden.txt\"\n"\
            "\n"\
            "## The maximum number of runs (invocations of targetRunner) that will performed. It\n"\
            "## determines the (maximum) budget of experiments for the tuning.\n"\
            "maxExperiments = ***\n"\
            "\n"\
            "## Indicates the number of decimal places to be considered for the\n"\
            "## real parameters.\n"\
            "digits = 2\n"\
            "\n"\
            "## A value of 0 silences all debug messages. Higher values provide\n"\
            "## more verbose debug messages.\n"\
            "debugLevel = 1\n"\
            "\n"\
            "## File to save tuning results as an R dataset, either absolute path\n"\
            "## or relative to execDir.\n"\
            "logFile = ***\n"\
            "\n"\
            "## Seed of the random number generator (must be a positive integer, NA\n"\
            "## means use a random seed).\n"\
            "seed = 11235813\n"\
            "\n"\
            "## Number of calls to targetRunner to execute in parallel. Less than 2\n"\
            "## means calls to targetRunner are sequentially executed.\n"\
            "parallel = ***\n"\
            "\n"\
            "## END of scenario file\n"\
            "############################################################################\n"

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("tsp_base", choices=["eil51", "pr107", "a280", "dsj1000"])
    parser.add_argument("item_relation_type", choices=["bsc", "unc", "usw"])
    parser.add_argument("number_of_items_per_city", choices=["1", "3", "5", "10"])
    parser.add_argument("max_experiments", help = "(maximum number of runs that are performed)")
    args = parser.parse_args()

    if os.path.isdir("train_instances"):
        if len(os.listdir("train_instances/")) > 0:
            os.system("rm train_instances/*")    
    else:
        os.mkdir("train_instances")

    for r, d, f in os.walk("../../../instances"):
        for _f in f:
            if "%s_%02d_%s" % (args.tsp_base, int(args.number_of_items_per_city), args.item_relation_type) in _f:            
                os.system("cp %s train_instances/." % (os.path.join(r, _f)))

    f = open("scenario.txt", "w")
    f.write(scenario.replace("maxExperiments = ***", "maxExperiments = %s" % (args.max_experiments))
                    .replace("logFile = ***", "logFile = \"log_%s_%02d_%s.Rdata\"" % (args.tsp_base, int(args.number_of_items_per_city), args.item_relation_type))
                    #.replace("parallel = ***", "parallel = %d" % (os.cpu_count()))
                    .replace("parallel = ***", "parallel = 16")
           )
    f.close()

    os.system("make -C ../")

    os.system("$IRACE_HOME/bin/irace > log_%s_%02d_%s.txt" % (args.tsp_base, int(args.number_of_items_per_city), args.item_relation_type))

    os.system("rm train_instances/*")
    os.system("scenario.txt")

