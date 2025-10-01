#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import math

EPS = 10e-4

def check_solution(input_file, solution_file, draw=False):

    file = open(input_file, 'r')
    lines = file.readlines()

    # lines[0]  // PROBLEM NAME:    a280-ThOP
    # lines[1]  // KNAPSACK DATA TYPE: bounded strongly corr
    # lines[2]  // DIMENSION:  280
    number_of_cities = int(lines[2].split()[-1])
    # lines[3]  // NUMBER OF ITEMS: 279
    number_of_items = int(lines[3].split()[-1])
    # lines[4]  // CAPACITY  OF  KNAPSACK: 25936
    capacity_of_knapsack = int(lines[4].split()[-1])
    # lines[5]  // MAX TIME: xxx
    max_time = float(lines[5].split()[-1])
    # lines[6]  // MIN SPEED: 0.1
    min_speed = float(lines[6].split()[-1])
    # lines[7]  // MAX SPEED: 1
    max_speed = float(lines[7].split()[-1])
    # lines[8]  // EDGE_WEIGHT_TYPE: CEIL_2D
    # lines[9]  // NODE_COORD_SECTION(INDEX, X, Y):

    vertices = {}
    for line in lines[10:10+number_of_cities]:
        index, pos_x, pos_y = line.split()
        vertices[int(index)] = (float(pos_x), float(pos_y), )

    # lines[10+number_of_cities]  // ITEMS SECTION(INDEX, PROFIT, WEIGHT, ASSIGNED NODE NUMBER):

    items = {}
    for line in lines[11+number_of_cities:11+number_of_cities+number_of_items]:
        index, profit, weight, id_city = line.split()
        items[int(index)] = (float(profit), float(weight), int(id_city), )

    file.close()

    is_feasible = True
    error_info = ""

    file = open(solution_file, 'r')
    lines = file.readlines()
    lines[0] = lines[0].replace(',',' ').replace('[','').replace(']','').split()
    performed_tour = [1] + [int(city) for city in lines[0]] + [number_of_cities]

    for city in performed_tour:
        if city < 1 or city > number_of_cities:
            is_feasible = False
            error_info += "[some visited city does not exist] "
            break

    lines[1] = lines[1].replace(',',' ').replace('[','').replace(']','').split()
    collected_items = [int(item) for item in lines[1]]

    for item in collected_items:
        if item < 1 or item > number_of_items:
            is_feasible = False
            error_info += "[some collected item does not exist] "
            break

    set_items = {item:1 for item in collected_items}
    if len(set_items) != len(collected_items):
        is_feasible = False
        error_info += "[there are repeated items] "

    for item in collected_items:
        if 1 <= item <= number_of_items:
            was_collected = False
            for city in performed_tour:
                if 1 <= city <= number_of_cities:
                    if items[item][2] == city:
                        was_collected = True
                        break
            if not was_collected:
                is_feasible = False
                error_info += "[the city of some collected item is not visited] "
                break


    if not is_feasible:
        print("%-30s %-30s %s" %(input_file, solution_file, "infeasible solution " + str(error_info)))
        return

    accumulate_profit = [0 for _ in range(0, len(performed_tour) + 1)]
    accumulate_weight = [0 for _ in range(0, len(performed_tour) + 1)]
    for i in range(len(performed_tour)):
        for j in range(len(collected_items)):
            if items[collected_items[j]][2] == performed_tour[i]:
                accumulate_profit[i] += items[collected_items[j]][0]
                accumulate_weight[i] += items[collected_items[j]][1]

    v = (max_speed - min_speed) / capacity_of_knapsack

    prev = 1
    current_capacity_of_knapsack = 0
    current_time = 0.0
    current_profit = 0.0
    for i in range(1, len(performed_tour)):
        current_time += math.ceil(math.sqrt(sum((vertices[prev][k]-vertices[performed_tour[i]][k])**2 for k in range(2)))) / (max_speed - v * current_capacity_of_knapsack)
        # print("%d %d %.5f %.5f" % (prev, performed_tour[i], current_time, max_time))
        current_capacity_of_knapsack += accumulate_weight[i]
        current_profit += accumulate_profit[i]
        prev = performed_tour[i]        

    if current_capacity_of_knapsack - EPS > capacity_of_knapsack:
        is_feasible = False
        error_info += "[capacity exceeded] "

    if current_time - EPS > max_time:
        is_feasible = False
        error_info += "[time exceeded] %f %f " % (current_time, max_time)

    if is_feasible:
        print("%-30s %-30s %.3f" % (input_file, solution_file, current_profit))
    else:
        print("%-30s %-30s %s" %(input_file, solution_file, "infeasible solution " + str(error_info)))
        return

    if not draw:
        return

    import matplotlib.pyplot as plt

    x_min = y_min = 10000000
    x_max = y_max = -1

    for key, value in vertices.items():
        x_min = min(x_min, value[0])
        x_max = max(x_max, value[0])
        y_min = min(y_min, value[1])
        y_max = max(y_max, value[1])
    x_min = int(math.floor(x_min - 0.10 * math.fabs(x_max)))
    x_max = int(math.ceil(x_max + 0.10 * math.fabs(x_max)))
    y_min = int(math.floor(y_min - 0.10 * math.fabs(y_max)))
    y_max = int(math.ceil(y_max + 0.10 * math.fabs(y_max)))

    plt.xlim(x_min, x_max), plt.ylim(y_min, y_max)

    x = [vertices[i][0] for i in performed_tour]
    y = [vertices[i][1] for i in performed_tour]

    current_capacity_of_knapsack = 0

    width = [1]
    for i in range(1, len(performed_tour)):
        current_capacity_of_knapsack += accumulate_weight[i]
        width.append(1 + 4 * (current_capacity_of_knapsack / capacity_of_knapsack))

    for i in range(len(x) - 1):
        plt.plot(x[i:i + 2], y[i:i + 2], linewidth=width[i], color = 'black')

    x_all = [value[0] for key, value in vertices.items()]
    y_all = [value[1] for key, value in vertices.items()]
    point_size_all = []
    for key, value in vertices.items():
        _sum_profit = 0
        for item_key, item_value in items.items():
            if item_value[2] == key:
                _sum_profit += item_value[0]
        point_size_all.append(_sum_profit)

    point_max_profit = max(point_size_all)

    for i in range(len(point_size_all)):
        point_size_all[i] = 10.0 + 50.0 * (point_size_all[i] / point_max_profit)

    print(point_size_all)

    plt.scatter(x_all, y_all, marker = 'o', s=point_size_all, color = 'black')

    source = [vertices[1]]
    plt.plot(*zip(*source), marker = 'o', markersize = 8.0, color = 'white')
    plt.plot(*zip(*source), marker = '>', markersize = 8.0, color = 'green')

    destination = [vertices[number_of_cities]]
    plt.plot(*zip(*destination), marker = 'o', markersize = 7.0, color = 'white')
    plt.plot(*zip(*destination), marker = 's', markersize = 7.0, color = 'red')

    plt.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    plt.xticks([])
    plt.yticks([])

    plt.savefig(solution_file + ".pdf", bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="")
    parser.add_argument("solution_file", help="")
    parser.add_argument("-d", "--draw", help="save a graphical representation of the (feasible) solution [matplotlib.pyplot is required]", action="store_true")
    args = parser.parse_args()
    check_solution(args.input_file, args.solution_file, args.draw)
