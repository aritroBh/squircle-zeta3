import os
base = os.path.dirname(os.getcwd())
d8 = os.path.join(base,'OCT8_symbol_test')
with open(os.path.join(d8,'step_14.py')) as f:
    print(f.read())
