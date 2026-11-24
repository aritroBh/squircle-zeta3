base = "<project>/research/"
for fn in ["planner/planner_check_naive.py", "planner/planner_check_F4.py"]:
    p = base + "08_theory/COR5_corner_classicality/" + fn
    with open(p) as f:
        print("="*80)
        print(fn)
        print("="*80)
        print(f.read())
