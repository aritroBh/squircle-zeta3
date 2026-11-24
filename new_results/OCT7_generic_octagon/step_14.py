import os
base="research/08_theory"
o6=os.path.join(base,"OCT6_octagon_closed")
print("oct5 subdir:", os.listdir(os.path.join(o6,"oct5")))
for s in ["step_08.py","step_09.py"]:
    print("="*30, s, "="*30)
    print(open(os.path.join(o6,s)).read())
