# A4: (1) L_simplified cache == tow2lib.L(r) (the PROVED construction), r=1..8, symbolic in (d,m) over Q(i) (sqrt2 pieces reduced),
#         plus a high-precision numeric cross-check at random (d,m);  (2) tow2lib.L(r) == my first-principles L_mine.
import sys, time, pickle, random
sys.dont_write_bytecode=True
W='<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, W)
import sympy as sp, mpmath as mp
import tow2lib as T2
d,m = T2.d,T2.m; I=sp.I
cache = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open(W+'/L_simplified.pkl','rb')).items()}
Lm = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('L_mine.pkl','rb')).items()}
mp.mp.dps = 50
random.seed(3)
t0=time.time()
for r in range(1,9):
    Lr = T2.L(r)
    has_sqrt = Lr.has(sp.sqrt(2))
    diff = sp.cancel(sp.together(Lr - cache[r]))
    diff2 = sp.cancel(sp.together(Lr - Lm[r]))
    # numeric check at 3 random points (d,m reals), 50 digits
    mx = 0
    for _ in range(3):
        dv = sp.Rational(random.randint(10,60),1)+sp.Rational(random.randint(1,9),10); mv = sp.Rational(random.randint(1,30),1)+sp.Rational(1,7)
        a = sp.N(Lr.subs({d:dv,m:mv}), 50); b = sp.N(cache[r].subs({d:dv,m:mv}), 50)
        mx = max(mx, abs(complex(a-b)))
    print(f'r={r}: T2.L(r) contains sqrt2 before reduction: {has_sqrt}; T2.L(r) - cache = {diff}; T2.L(r) - L_mine = {diff2}; numeric max|diff| (3 pts) = {mx:.1e}  ({time.time()-t0:.0f}s)', flush=True)
