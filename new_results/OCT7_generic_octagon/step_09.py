import os
base="research/08_theory"
o5=os.path.join(base,"OCT5_octagon_exact")
for s in ["step_01.py","step_02.py","step_03.py"]:
    print("="*30, s, "="*30)
    print(open(os.path.join(o5,s)).read())
