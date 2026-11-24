import os
M4E='<project>/research/08_theory/M4e_w5_lchi'
for f in ['w5_symbol.py','w5_search.py','w5_cobracket.py','lchi_valrel.py','nice_sat_s0.py','display_nice.py','planner/planner_reduce.py']:
    print(f, os.path.getsize(f'{M4E}/{f}'))
print()
print(open(f'{M4E}/planner/planner_reduce.py').read())
