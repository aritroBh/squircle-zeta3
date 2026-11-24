import os
base = '../PGEN1R_referee_alphabet'
for f in ['step_70.py','step_71.py','step_72.py','step_73.py','step_74.py','step_75.py','step_76.py']:
    print(f"===== {f} =====")
    print(open(os.path.join(base,f)).read())
