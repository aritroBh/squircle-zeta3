import os
base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon')
with open(os.path.join(d1,'step_60.py')) as f:
    print(f.read())
