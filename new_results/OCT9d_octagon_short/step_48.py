import sys, time, pickle
sys.path.insert(0, ".")
from oct9d_common import *
TW = TW3()
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
sel = ['-9/16', '4/9', '1/5', '1/3', '-5/3', '-1/5', '5/4', '2/5']
selq = [sp.Rational(g) for g in sel]
PB = [(sp.Rational(w),p) for w,p in [('-1/3',2),('-1/3',3),('-1/3',5),('-1/4',2),('-1/4',3),('-1/4',5),('1/4',2),('1/4',3),('1/4',5),
      ('-5/3',2),('-5/3',3),('-5/3',5),('-5/4',2),('-5/4',3),('-5/4',5),('1/5',3),('1/5',5),('5/8',2),('5/8',3)]]
PMAT = sp.Matrix.hstack(*[col_prod(w,p,SS235) for w,p in PB])
M = sp.Matrix.hstack(*[Gsym[g] for g in selq], PMAT)
print("det [G8|PB] =", M.det(), "=", sp.factorint(abs(M.det())))
for e in (1,2,3):
    x = M.LUsolve(sp.Matrix(TW[e]))
    assert M*x == sp.Matrix(TW[e])
    cl = [sp.Rational(c) for c in x[:8]]; cp = [sp.Rational(c) for c in x[8:]]
    print(f"--- e{e} ---")
    print(" Li3:", [str(c) for c in cl])
    print(" prod:", [str(c) for c in cp])
