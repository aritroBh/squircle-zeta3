import os
base="research/08_theory"
o5=os.path.join(base,"OCT5_octagon_exact")
for s in ["step_04.py","step_05.py","step_06.py","step_07.py","step_08.py"]:
    print("="*30, s, "="*30)
    print(open(os.path.join(o5,s)).read()[:3500])
    print()
