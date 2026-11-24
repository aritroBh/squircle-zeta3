import os, pickle, time, sys
import sympy as sp
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
d, m = sp.symbols('d m')
Ls_raw = pickle.load(open(work + '/L_simplified.pkl','rb'))
L = {k: sp.sympify(v, locals={'d':d,'m':m,'I':sp.I}) for k,v in Ls_raw.items()}

from collections import Counter
def parts(rem, minpart, cur):
    if rem == 0:
        yield tuple(sorted(cur)); return
    for p in range(minpart, rem+1):
        yield from parts(rem-p, p, cur+[p])

def bell(delta):
    tot = sp.Integer(0)
    for part in parts(delta,1,[]):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz,mult in C.items():
            coef *= (L[sz]**mult)/sp.factorial(mult)
        tot += coef
    return tot

for delta in [7,8]:
    t0 = time.time()
    a = bell(delta)
    a = sp.cancel(a)
    num = sp.numer(a); den = sp.denom(a)
    print(f"a_{delta}: numdeg_d={sp.Poly(num,d).degree()} terms={len(sp.Add.make_args(sp.expand(num)))} denom={sp.factor(den)}  ({round(time.time()-t0,1)}s)", flush=True)
    pickle.dump(str(a), open(os.path.join(work, f'a{delta}_closed.pkl'),'wb'))
print("done")
