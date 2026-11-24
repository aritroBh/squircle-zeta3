import os
prev = "<project>/research/08_theory/COR3_classical_proof/"
if os.path.exists(prev):
    print(sorted(os.listdir(prev)))
else:
    print("not found")
