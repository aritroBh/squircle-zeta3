# A2: symbolic, in (d,m).  Independent first-principles a_delta (series of A(i+x)/A(i), h = x^2/(4(i+x))), delta=1..8:
#   a_delta^{mine} = [x^delta] (1+x/i)^{2d-1} (1+x/2i)^{2-2d} (8/G(i+x)) (1 + R1 h + R2 h^2/2 + R3 h^3/6 + R4 h^4/24),  R_l = N_l/Pi_l  (checked vs S definition in a1_R.py)
#   compare with   (1) tow4lib.arho(delta) = Bell(L cache),  (2) a7_closed.pkl / a8_closed.pkl,  (3) the displayed N_7 in proof_depth8.md
#   and L_r^{mine} (log of the series) == L_simplified cache == tow2lib.L(r), r=1..8.
import sys, time, pickle, re
sys.dont_write_bytecode=True
W='<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, W)
import sympy as sp
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
I = sp.I
x = sp.symbols('x')
r1,r2,r3,r4 = sp.symbols('r1:5')
Nn = 8
def trunc(e):
    e = sp.expand(e)
    return sum(e.coeff(x,k)*x**k for k in range(Nn+1))
def sbinom_series(c, alpha):      # (1 + c x)^alpha, alpha polynomial in d
    out=0; cp=1; coef=sp.Integer(1)
    for n in range(Nn+1):
        out += sp.expand(coef*cp)*x**n
        coef = sp.expand(coef*(alpha-n)/(n+1)); cp = cp*c
    return out
f1 = sbinom_series(-I, 2*d-1)
f2 = sbinom_series(-I/2, 2-2*d)
v = I + x
Gp = sp.expand(v**4 - 6*v**2 + 1)
assert Gp.coeff(x,0) == 8
Gser = sp.expand(Gp/8)
# inverse series of Gser
inv = [sp.Integer(1)]
gc = [sp.expand(Gser).coeff(x,k) for k in range(Nn+1)]
for n in range(1,Nn+1):
    inv.append(sp.expand(-sum(gc[k]*inv[n-k] for k in range(1,n+1))))
Ginv = sum(inv[k]*x**k for k in range(Nn+1))
# h = x^2/(4(i+x)) = (x^2/(4i)) * (1 + x/i)^{-1}
hser = trunc(x**2/(4*I) * sum((I*x)**k for k in range(Nn+1)))   # (1+x/i)^{-1} = sum (i x)^k  [BUG FIX: was (-i x)^k]
assert sp.simplify(sp.series(hser - x**2/(4*(I+x)), x, 0, Nn+1).removeO()) == 0
# check h = t(i+x)-t* directly
vv = I + x
tchk = sp.series(sp.simplify((vv**2+2*vv-1)/(4*vv) - (1+I)/2 - x**2/(4*(I+x))), x, 0, 10).removeO()
print('h(x) = t(i+x) - t* identity:', sp.simplify(tchk) == 0)
hp = [sp.Integer(1)]
for l in range(1,5): hp.append(trunc(hp[-1]*hser))
Sr = 1 + r1*hp[1] + r2*hp[2]/2 + r3*hp[3]/6 + r4*hp[4]/24
tot = trunc(trunc(trunc(f1*f2)*Ginv)*Sr)
asym = {k: sp.expand(tot.coeff(x,k)) for k in range(Nn+1)}
subsR = {r1: T2.R[1], r2: T2.R[2], r3: T2.R[3], r4: T2.R[4]}
mine = {0: sp.Integer(1)}
t0=time.time()
for k in range(1,Nn+1):
    mine[k] = sp.cancel(sp.together(asym[k].subs(subsR, simultaneous=True)))
    print(f'a_{k}^mine built ({time.time()-t0:.1f}s); r-degree structure: linear in R_l (max R-power {sp.Poly(asym[k], r1,r2,r3,r4).total_degree()})', flush=True)
# (1) vs tow4lib.arho (Bell of L cache)
res = {}
for k in range(1,Nn+1):
    res[k] = sp.cancel(sp.together(mine[k] - T4.arho(k)))
    print(f'a_{k}: mine - tow4lib.arho (Bell of L cache) = {res[k]}   ({time.time()-t0:.1f}s)', flush=True)
# (2) pickles
for k,fn in ((7,'a7_closed.pkl'),(8,'a8_closed.pkl')):
    raw = pickle.load(open(W+'/'+fn,'rb'))
    obj = raw if not isinstance(raw,dict) else raw
    print(fn, type(raw), (list(raw.keys()) if isinstance(raw,dict) else ''))
pickle.dump({k:str(mine[k]) for k in mine}, open('a_mine.pkl','wb'))
# L_r^mine via inverse Bell recursion
Lm = {}
for n in range(1,Nn+1):
    Lm[n] = sp.cancel(sp.together(mine[n] - sp.Rational(1,n)*sum(r*Lm[r]*mine[n-r] for r in range(1,n))))
cache = {k: sp.sympify(vv_, locals={'d':d,'m':m,'I':I}) for k,vv_ in pickle.load(open(W+'/L_simplified.pkl','rb')).items()}
for n in range(1,Nn+1):
    dd = sp.cancel(sp.together(Lm[n]-cache[n]))
    print(f'L_{n}: mine(log of first-principles series) - L_simplified cache = {dd}', flush=True)
pickle.dump({k:str(Lm[k]) for k in Lm}, open('L_mine.pkl','wb'))
print('done', time.time()-t0)
