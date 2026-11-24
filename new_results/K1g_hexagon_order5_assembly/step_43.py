import os
base = '<project>/research/08_theory'
k1e = os.path.join(base, 'K1e_hexagon_closed')
p = os.path.join(k1e, 'planner_check_irrational.py')
if os.path.exists(p):
    with open(p) as f: print(f.read())
else:
    print("not found in K1e; searching...")
    for d in os.listdir(base):
        cand = os.path.join(base, d, 'planner_check_irrational.py')
        if os.path.exists(cand):
            print("FOUND in", d)
