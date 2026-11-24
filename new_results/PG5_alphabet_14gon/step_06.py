import os
base = '../PGEN1R_referee_alphabet'
for f in ['step_67.py','step_68.py']:
    print(f"===== {f} =====")
    print(open(os.path.join(base,f)).read())
