base = '<project>/research/08_theory/'
d3 = base+'TOW3_tower_depth6/'
for s in ['step_28.py','step_29.py','step_30.py','step_31.py','step_32.py']:
    t = open(d3+s).read()
    print('='*20, s, len(t))
    print(t[:2500])
