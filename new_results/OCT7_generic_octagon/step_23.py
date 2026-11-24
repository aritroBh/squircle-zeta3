import os
base="research/08_theory"
o5=os.path.join(base,"OCT5_octagon_exact")
print(open(os.path.join(o5,"inputs","OCT4_result.md")).read()[:7000])
