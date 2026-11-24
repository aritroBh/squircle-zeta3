import sys, json
from id_lib import *
dps = 300; set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
I = load_I(dps)
c = [combo(I, f"c{j}") for j in range(3)]
A = atoms(); D = 285
def show(name, tgt, names, D=D):
    S = with_sqrt2({n: A[n] for n in names})
    vals = [tgt] + list(S.values()); labels = ["target"] + list(S.keys())
    res = find_relation(vals, D)
    co, s = res[0]
    print(f"{name}: n={len(vals)} D={D} best row max|coef|={max(abs(x) for x in co)}  residual={mp.nstr(s,4)}  target coef={co[0]}  mp.dps={mp.mp.dps}")
    return co, labels
mz = ["lam^4", "pi2 lam^2", "pi4", "z3 lam", "Li4(1/2)"]
for nm, t in (("c2", c[2]), ("c1", c[1]), ("c0", c[0])):
    co, lab = show(nm + " vs alt-MZV4 x{1,s2}", t, mz)
