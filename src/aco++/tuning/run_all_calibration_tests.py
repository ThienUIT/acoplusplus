import os

tsp_base = ["eil51", "pr107", "a280", "dsj1000"]
item_relation_type = ["bsc", "unc", "usw"]
number_of_items_per_city = [1, 3, 5, 10]

for tsp in tsp_base:
    for num in number_of_items_per_city:
        for rel in item_relation_type:
            os.system("python run_irace.py %s %s %d 5000" % (tsp, rel, num))

