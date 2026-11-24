import sys
sys.path.insert(0, ".")
from qs3_symbol import *
print("constants and their logvecs:")
for nm, z in [("alpha=beta=(2-r3)/2", q3(1, F(-1,2))), ("c0_0=2-r3", q3(2,-1)),
              ("c0_1=(2r3-3)/2", q3(F(-3,2),1)), ("c0_2=(2-r3)/2", q3(1,F(-1,2))),
              ("kappa0", (F(-1),F(1,2))), ("kappa1", (F(-1),F(2,3))), ("kappa2",(F(-1),F(1))), ("kappa3",(F(-1),F(1,3))),
              ("2", q3(2)), ("3", q3(3)), ("1/2", q3(F(1,2))), ("1/3", q3(F(1,3))),
              ("L/u=2-r3", q3(2,-1)), ("s0/u", q3(-1,F(1,2))), ("u^2=2-r3", q3(2,-1))]:
    print(f"  {nm}: {logvec(z)}")
# check the itab progress
print("\n--- itab log ---")
print(open("itab_build.log").read())
