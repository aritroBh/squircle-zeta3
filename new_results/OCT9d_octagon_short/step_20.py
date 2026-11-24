import sys, time
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
# Build [G|D] matrix for edges 1-3 and verify stored sols (PROVED)
cols = [col_Li3(g, SS235) for g in Gargs] + [col_prod(w,p,SS235) for (w,p) in Dargs]
M235 = sp.Matrix.hstack(*cols)
print("M235", M235.rows, "x", M235.cols, f"({time.time()-t0:.1f}s)")
rk = M235.rank(); print("rank [G|D] =", rk, f"({time.time()-t0:.1f}s)")
G_only = sp.Matrix.hstack(*cols[:99]); print("rank G (99 Li3 alone) =", G_only.rank())
for e in (1,2,3):
    sol = sp.Matrix([sp.Rational(s) for s in D['s123'][e]['sol']])
    print(f"edge {e}: [G|D]*sol == TW3 ?", M235*sol == sp.Matrix(TW[e]))
    print(f"        TW3 in span G alone? rank G == rank [G|t]:", G_only.rank(), (sp.Matrix.hstack(G_only, sp.Matrix(TW[e]))).rank())
