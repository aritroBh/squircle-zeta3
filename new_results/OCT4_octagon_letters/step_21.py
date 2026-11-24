import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp
from itertools import permutations
sqrt2 = sp.sqrt(2); I = sp.I
oo = sp.oo
# letter set
tstar = (1+I*(1+sqrt2))/2
letters = {'0':sp.Integer(0),'1':sp.Integer(1),'inf':oo,
           '-(1+r2)':-(1+sqrt2),'2+r2':2+sqrt2,'t*':tstar,'conj':sp.conjugate(tstar)}
def cr(a,b,c,d):
    # (a-c)(b-d)/((a-d)(b-c)) with infinity handling via limit
    def sub(expr):
        return expr
    def term(x,y):
        if x is oo and y is oo: return sp.Integer(1)
        if x is oo: return sp.Integer(1)  # (oo - y) ~ oo, handled by ratio
        if y is oo: return sp.Integer(-1)
        return x-y
    # use formula with symbolic limit: cr = ((a-c)/(a-d))/((b-c)/(b-d))
    def ratio(p,q,r,ss):
        # (p-q)/(r-ss)
        num = term(p,q); den = term(r,ss)
        if num==sp.Integer(1) and den==sp.Integer(1): return sp.Integer(1)
        return num/den
    return sp.simplify(ratio(a,c,a,d)/ratio(b,c,b,d))

# Real-letter cross-ratios first
real_keys = ['0','1','inf','-(1+r2)','2+r2']
realset=set()
for combo in permutations(real_keys,4):
    v = cr(*(letters[k] for k in combo))
    realset.add(sp.simplify(v))
print("Distinct real-letter cross-ratios:", len(realset))
for v in sorted(realset, key=lambda z: float(sp.N(sp.re(z))) if not z.has(I) else 999):
    print("   ", sp.nsimplify(v), "=", sp.N(v,8))
log("step: real-letter cross-ratios computed, %d distinct"%len(realset))
