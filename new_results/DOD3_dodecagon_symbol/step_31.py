import sys, os
base = os.path.dirname(os.getcwd())
d2 = os.path.join(base,'DOD2_dodecagon_finish')
with open(os.path.join(d2,'step_14.py')) as f:
    print(f.read())
