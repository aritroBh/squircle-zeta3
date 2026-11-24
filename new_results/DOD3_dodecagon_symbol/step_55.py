import os
base = os.path.dirname(os.getcwd())
d8 = os.path.join(base,'OCT8_symbol_test')
for f in ['step_32.py','step_43.py']:
    print('='*30, f)
    with open(os.path.join(d8,f)) as fh:
        print(fh.read()[:6000])
