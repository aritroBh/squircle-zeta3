import os
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
print(os.listdir(O9B))
src = open(os.path.join(O9B,"syminfra.py")).read()
print(src[:9000])
