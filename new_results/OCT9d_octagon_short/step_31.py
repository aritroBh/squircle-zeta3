import sys, pickle
sys.path.insert(0, ".")
from oct9d_common import *
D = load_all()
for e in (1,2,3):
    sol = [sp.Rational(s) for s in D['s123'][e]['sol']]
    nzG = [(Gargs[i], sol[i]) for i in range(99) if sol[i] != 0]
    nzD = [(Dargs[i], sol[99+i]) for i in range(27) if sol[99+i] != 0]
    print(f"edge {e}: nonzero Li3 = {len(nzG)}, nonzero prod = {len(nzD)}")
    print("  Li3 args used:", [str(g) for g,c in nzG])
    print("  prod used:", [(str(w),p) for (w,p),c in nzD])
sol0 = [sp.Rational(s) for s in D['s0']['solB0']]
n7 = len(sols7); nw = len(negw7)*4; nl = len(logprods7)
nzG0 = [(sols7[i], sol0[i]) for i in range(n7) if sol0[i] != 0]
nzD0 = [((negw7[i//4], (2,3,5,7)[i%4]), sol0[n7+i]) for i in range(nw) if sol0[n7+i] != 0]
nzL0 = [(logprods7[i], sol0[n7+nw+i]) for i in range(nl) if sol0[n7+nw+i] != 0]
print(f"edge 0: nonzero Li3 = {len(nzG0)}, prod = {len(nzD0)}, log3 = {len(nzL0)}")
print("  Li3 args:", [str(g) for g,c in nzG0])
print("  prods:", [(str(w),p) for (w,p),c in nzD0], " log3:", [t for t,c in nzL0])
