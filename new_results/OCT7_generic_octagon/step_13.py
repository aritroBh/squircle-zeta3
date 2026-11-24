import os
base="research/08_theory"
o6=os.path.join(base,"OCT6_octagon_closed")
# OCT6 step_01 - monomials; step_03 - exact IBP engine
for s in ["step_01.py"]:
    print("="*30, s, "="*30)
    print(open(os.path.join(o6,s)).read())
