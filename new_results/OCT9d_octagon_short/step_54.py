import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
import flint
t0=time.time()
TW = TW3()

# ---- final bases ----
G8 = [sp.Rational(g) for g in ['-9/16','4/9','1/5','1/3','-5/3','-1/5','5/4','2/5']]
PB19 = [(sp.Rational(w),p) for w,p in [('-1/3',2),('-1/3',3),('-1/3',5),('-1/4',2),('-1/4',3),('-1/4',5),
        ('1/4',2),('1/4',3),('1/4',5),('-5/3',2),('-5/3',3),('-5/3',5),('-5/4',2),('-5/4',3),('-5/4',5),
        ('1/5',3),('1/5',5),('5/8',2),('5/8',3)]]
B0 = [sp.Rational(g) for g in pickle.load(open("basis0.pkl","rb"))['basis0']]
PB0 = [(sp.Rational(w),p) for w,p in pickle.load(open("PB0.pkl","rb"))['PB0']]

def to_fmpq(cols, t):
    n = cols[0].rows
    M = flint.fmpq_mat(n, len(cols))
    for j, c in enumerate(cols):
        for i in range(n):
            x = sp.Rational(c[i]); M[i,j] = flint.fmpq(x.p, x.q)
    tf = flint.fmpq_mat(n,1)
    for i in range(n):
        x = sp.Rational(t[i]); tf[i,0] = flint.fmpq(x.p, x.q)
    return M, tf

# {2,3,5} solutions
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in G8}
PBcols = [col_prod(w,p,SS235) for w,p in PB19]
sols123 = {}
for e in (1,2,3):
    Mf, tf = to_fmpq([Gsym[g] for g in G8] + PBcols, TW[e])
    assert Mf.det() != 0
    xf = Mf.solve(tf)
    assert Mf*xf == tf
    sols123[e] = [sp.Rational(sp.Integer(int(xf[i,0].numerator)), sp.Integer(int(xf[i,0].denominator))) for i in range(27)]
    print(f"edge {e}: exact solve verified (PROVED), det = {Mf.det()}  ({time.time()-t0:.1f}s)")
# edge 0
Gsym0 = {g: col_Li3(sp.Rational(g), SS2357) for g in B0}
PB0cols = [col_prod(w,p,SS2357) for w,p in PB0]
Mf, tf = to_fmpq([Gsym0[g] for g in B0] + PB0cols, TW[0])
assert Mf.det() != 0
xf = Mf.solve(tf)
assert Mf*xf == tf
sol0 = [sp.Rational(sp.Integer(int(xf[i,0].numerator)), sp.Integer(int(xf[i,0].denominator))) for i in range(64)]
print(f"edge 0: exact solve verified (PROVED), det = {Mf.det()}  ({time.time()-t0:.1f}s)")
pickle.dump({'G8': [str(g) for g in G8], 'PB19': [(str(w),p) for w,p in PB19],
             'B0': [str(g) for g in B0], 'PB0': [(str(w),p) for w,p in PB0],
             'sols123': {e: [str(c) for c in sols123[e]] for e in sols123},
             'sol0': [str(c) for c in sol0]}, open("short_form.pkl","wb"))
print("saved short_form.pkl")
