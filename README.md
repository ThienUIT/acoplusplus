# Efficiently solving the thief orienteering problem with a max-min ant colony optimization algorithm

This repository contains the source code and data associated to the paper ["Efficiently solving the thief orienteering problem with a max-min ant colony optimization algorithm"](https://link.springer.com/article/10.1007/s11590-021-01824-y) by Jonatas B. C. Chagas and Markus Wagner. The paper presents a swarm intelligence based on ant colony optimization combined with a packing heuristic algorithm for solving the Thief Orienteering Problem (ThOP). Our algorithm (called ACO++) is an extension of our previously proposed [ACO algorithm](https://www.sciencedirect.com/science/article/abs/pii/S0167637720301255). In addition to the ACO++ code, we make available in this repository the codes of 3 algorithms proposed in the literature for the ThOP: [ILS and BRKGA](https://ieeexplore.ieee.org/document/8477853), and [ACO](https://www.sciencedirect.com/science/article/abs/pii/S0167637720301255).

### Compiling the code

Before running our ACO++ algorithm, it is needed to compile its code. To this end, just run the following command:

```console
$ make
```

### Usage:

```console
$ ./acothop [parameters]

Parameters:

  -i, --inputfile       inputfile (ThOP format necessary)
  -o, --outputfile      outputfile
  -m, --ants            number of ants
  -a, --alpha           influence of pheromone trails
  -b, --beta            influence of heuristic information
  -e, --rho             pheromone trail evaporation
  -p, --ptries          number of tries to construct a packing plan from a give tour
  -l, --localsearch     0: no local search   1: 2-opt   2: 2.5-opt   3: 3-opt
  -t, --time            maximum time for each trial  
      --seed            seed for the random number generator
      --log             save an extra file (<outputfile>.log) with log messages

```

We provide a python script (see "src/aco++/run_aco++_experiments.py") for running all the computational experiments reported in our paper concerning our ACO++. In addition, with the same purpose, there is a python script for running each of the other algorithms available here.
