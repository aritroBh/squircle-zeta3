base = "<project>/research/08_theory/"
import os
# find the evaluator
print(os.listdir(base+"M4e_w5_lchi/planner_main"))
with open(base+"M4e_w5_lchi/planner_main/planner_main_check_best.py") as f:
    print(f.read())
