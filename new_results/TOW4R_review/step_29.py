
import os, sys, pickle, time
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2
d, m = T2.d, T2.m
i = sp.I
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1
# my own independent computation of tau1(8), tau2(8) from the DL definition in exact Q(sqrt2, i)
def tau1_def(n):
    s = n+1
    return (-1)**s * ( (i-bq)**(-(s+1)) + (i+aq)**(-(s+1)) - i**(-(s+1)) )
def tau2_def(n):
    s = n+1
    return (-1)**s * ( (i-aq)**(-(s+1)) + (i+bq)**(-(s+1)) - i**(-(s+1)) )
t1c = sp.Rational(-169,128) + sp.Rational(41,128)*i
t2c = sp.Rational(-169,128) - sp.Rational(41,128)*i
d1 = sp.simplify(sp.expand(tau1_def(8) - t1c))
d2 = sp.simplify(sp.expand(tau2_def(8) - t2c))
print("tau1(8) exact diff:", d1, " tau2(8) exact diff:", d2)
# 40-digit numeric cross-check
from mpmath import mp
mp.dps = 40
ii = mpc = mp.mpc(0,1); r2 = mp.sqrt(2)
a = r2+1; b = r2-1
def t1mp(n):
    s = n+1
    return (-1)**s*( (ii-b)**(-(s+1)) + (ii+a)**(-(s+1)) - ii**(-(s+1)) )
print("tau1(8) mp:", mp.nstr(t1mp(8), 30))
print("claim    mp:", mp.nstr(mp.mpf(-169)/128 + mp.mpf(41)/128*ii, 30))
# conjugation symmetry at n=8 and full table n=0..8
sys.path.insert(0, work+"/tow4_copy")
import tow4lib as T4
ok = all(sp.simplify(T4.tau1f(n) - (-1)**n*sp.conjugate(T4.tau2f(n)))==0 for n in range(0,9))
print("conj(tau2_n)=(-1)^n tau1_n n=0..8:", ok)
ok2 = all(sp.simplify(T4.tau2f(n) - (-1)**n*sp.conjugate(T4.tau1f(n)))==0 for n in range(0,9))
print("conj(tau1_n)=(-1)^n tau2_n n=0..8:", ok2)
# towers pickle meta
H = pickle.load(open("<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl",'rb'))
print("_meta:", H.get('_meta'))
keys = [k for k in H if k != '_meta']
print("nkeys:", len(keys), "max degree:", max(sum(k) for k in keys))
import collections
degcount = collections.Counter(sum(k) for k in keys)
print("degrees present:", sorted(degcount.items())[:8], "...", sorted(degcount.items())[-4:])
# cell count for RHO delta=7,8
for delta in [7,8]:
    n = sum(1 for (p,q) in keys if 2*(p+q)-4 >= delta)
    print(f"keys with 2d-4>={delta}:", n)
