cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
import os
print(os.listdir(cor))
with open(cor+"step_01.py") as f:
    print(f.read()[:4000])
