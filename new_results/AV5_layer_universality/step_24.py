import sympy as sp, mpmath, pickle
from mpmath import mp, mpf, nstr
mp.dps = 90
mons = pickle.load(open('ckpt_monomials.pkl','rb'))
needI, needK2, needK3 = set(), set(), set()
def collect(a,b,c,d,e):
    if d==0 and e==0:
        if c==0: needI.add((a,b))
        elif c==1: needI.add((a+2,b+1))
        elif c==2:
            if b==1: needK2.add(a+6)
            elif b==2: needK3.add(a+6)
        return
    if c==0 and e==0 and d==1:
        if b==0: return
        needI.add((a+4,b+1)); needI.add((a+4,b+1))
        if b==1: pass
        else: needK2.add(a+6)
        return
    if c==1 and d==1 and e==0 and b==0: return
    if e==1 and b==0 and c==0 and d==0: return
    if e==1 and b==1 and c==0 and d==0:
        needI.add((a+6,2)); needI.add((a+6,2)); needI.add((a+6,2))
        return
    raise ValueError((a,b,c,d,e))
for k in range(1,9):
    for (cs,a,b,c,d,e) in mons[k]:
        collect(a,b,c,d,e)
print("need I:", sorted(needI))
print("need K2:", sorted(needK2))
print("need K3:", sorted(needK3))
