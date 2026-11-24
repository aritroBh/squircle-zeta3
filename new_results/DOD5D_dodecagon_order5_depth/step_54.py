import sys
sys.path.insert(0, ".")
from qs3_symbol import *
# the six G-letters: 0, 1, p0..p3
p0 = q3(4, 2)      # 4+2r3 = V^2
p1 = q3(-3, -2)    # -(3+2r3) = -r3*U
p2 = q3(F(-1,2), F(-1,2))   # -(1+r3)/2
p3 = q3(F(3,2), F(1,2))     # (3+r3)/2
LETS = [ZERO, ONE, p0, p1, p2, p3]
names = ["0", "1", "p0", "p1", "p2", "p3"]
print("letter factorisations log|p_i| = aX+bY+cZ+dW:")
for i, z in enumerate(LETS):
    if i < 2: continue
    print(f"  {names[i]} = {z}: {logvec(z)}")
print("\nall pairwise differences a-b (a row, b col):")
ok = True
for i, a in enumerate(LETS):
    for j, b in enumerate(LETS):
        if i == j: continue
        try:
            lv = logvec(q3sub(a, b))
            tag = "ok"
        except Exception as e:
            lv = None; tag = f"FAIL {e}"; ok = False
        print(f"  {names[i]}-{names[j]}: {lv}  {tag}")
print("all S-unit:", ok)
# constants appearing: alpha=(2-r3)/2, c0s, kappas
print("\nconstants:")
for nm, z in [("alpha=beta=(2-r3)/2", q3(1, F(-1,2))), ("c0_0=2-r3", q3(2,-1)),
              ("c0_1=(2r3-3)/2", q3(F(-3,2),1)), ("c0_2=(2-r3)/2", q3(1,F(-1,2))),
              ("kappa0", (F(-1),F(1,2))), ("kappa1", (F(-1),F(2,3))), ("kappa2",(F(-1),F(1))), ("kappa3",(F(-1),F(1,3))),
              ("2", q3(2)), ("3", q3(3)), ("L/u=2-r3", q3(2,-1)), ("h/u=1/2", q3(F(1,2))), ("1/3", q3(F(1,3)))]:
    print(f"  {nm}: {logvec(z)}")
