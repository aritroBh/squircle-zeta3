import os
WORK = os.getcwd()
for f in ['step_01.py','step_02.py','step_03.py']:
    print('='*30, f)
    with open(f) as fh:
        print(fh.read())
