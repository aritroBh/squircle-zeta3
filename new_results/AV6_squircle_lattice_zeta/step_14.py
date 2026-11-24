
from squircle import *
from mpmath import mp, zeta, dirichlet, mpf, nstr, fabs
import time
setdps(60)
def beta(x): return dirichlet(x, [0,1,0,-1])
def rel(a,b):
    d = fabs(a-b)
    return d/fabs(b) if d else mpf('1e-999')
rows = []
ps = [0.5, 0.75, 1.0, 1.05, 1.5, 2, 2.5, 3, 4, 6, 10, 25, 50]
ss = [3, 4, 6]
for p_ in ps:
    for s_ in ss:
        t0 = time.time()
        ZA_ = ZA(p_, s_)
        ZB_, iB = ZB(p_, s_, info=True)
        r = rel(ZA_, ZB_)
        note = ""
        if p_ == 1:
            r1 = rel(ZA_, 4*zeta(s_-1)); note = f"exact4z:{nstr(r1,3)}"
        if p_ == 2:
            r2 = rel(ZA_, 4*zeta(s_/2)*beta(s_/2)); note = f"exactEp:{nstr(r2,3)}"
        rows.append((p_, s_, ZA_, r, note, time.time()-t0, iB['delta'], iB['est1']))
        print(f"p={p_:<5} s={s_}  A-B={nstr(r,4):>10}  {note:>16} d={nstr(iB['delta'],3):>9} e1={nstr(iB['est1'],3):>9} {time.time()-t0:5.1f}s")
import pickle
with open("sweep.pkl","wb") as f:
    pickle.dump([(str(a),str(b),str(c)) for a,b,c,_,_,_,_,_ in rows], f)
