import os
base = '../PGEN1R_referee_alphabet'
for f in ['step_31.py','step_33.py','step_50.py','step_51.py','step_52.py','step_53.py','step_54.py','step_55.py']:
    print(f"===== {f} =====")
    print(open(os.path.join(base,f)).read())
