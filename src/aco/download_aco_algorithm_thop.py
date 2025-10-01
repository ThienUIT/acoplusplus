import os

os.system("git clone https://github.com/jonatasbcchagas/aco_thop.git")
os.system("rm aco_thop/src/run_all_experiments.py")
os.system("mv aco_thop/src/* .")
os.system("rm -rf santos_and_chagas_programs")
os.system("rm -rf aco_thop")
