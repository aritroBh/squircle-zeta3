"""Compare TOW4's E3 formula (copied scripts) at p=3 with the actual Laurent coefficients (orders 8..1, and the Taylor term 0)
of the direct and fold parts of Phi_Delta computed from exact data (towers_D8_K52 + my own PF arithmetic)."""
import sys, os, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, os.getcwd())
import sympy as sp
from chk_feed import *          # H, QSI, pf ops, gm_from_tower, hdict
import prove_E3b78 as PE
from fractions import Fraction as Fr
p = PE.p
# --- data side
def lagr_vec(xs, ys, x):
    s=[Fr(0)]*4
    for i,xi in enumerate(xs):
        L=Fr(1)
        for j,xj in enumerate(xs):
            if i!=j: L*=Fr(x-xj, xi-xj)
        for t in range(4): s[t]+=ys[i][t]*L
    return tuple(s)
Delta_raw = {}
for dl in range(9):
    P=dl//2; xs=list(range(7,7+P+1)); ys=[H[(3,3)][dl][k] for k in xs]
    st=lagr_vec(xs,ys,6)
    Delta_raw[dl]=tuple(a-b for a,b in zip(H[(3,3)][dl][6],st))
rho0 = QSI(*H[(3,3)][0][6])
hD = {10-dl: QSI(*Delta_raw[dl]) for dl in range(9)}
G = gm_from_tower(hD)
Pd = pf_scale(pf_mul_DL(G, DL2), Fr(-3))            # direct part
Pf = pf_scale(pf_fold(pf_mul_DL(G, DL1)), Fr(-3))   # fold part
def laurent(P, r):
    """coefficient of (v-i)^{-r}, r>=1 from ipos; r=0: Taylor constant of the regular part at v=i"""
    real, ipos, ineg, const = P
    if r>=1:
        return ipos[r-1] if r<=len(ipos) else QSI(0)
    from fastred import IQ
    tot = const
    for q,c in real.items(): tot = tot + c/(IQ - q)
    for j,f in enumerate(ineg, start=1): tot = tot + f/((IQ+IQ)**j)
    return tot
res=[]
for dl in range(0,9):
    rc = 8-dl
    res.append((dl, rc, laurent(Pd,rc)/rho0 if not rho0.is_zero() else None, laurent(Pf,rc)/rho0))
# --- formula side at p=3
terms_cache={}
out=[]
for dl in range(0,9):
    ts = PE.feed_terms(dl)
    D = sp.nsimplify(sum(ts[:dl+1]).subs(p,3)) if True else None
    F = sp.nsimplify(sum(ts[dl+1:]).subs(p,3))
    out.append((dl, sp.simplify(D), sp.simplify(F)))
def toS(q): return sp.Rational(q.p0.numerator,q.p0.denominator)+sp.sqrt(2)*sp.Rational(q.p1.numerator,q.p1.denominator)+sp.I*(sp.Rational(q.p2.numerator,q.p2.denominator)+sp.sqrt(2)*sp.Rational(q.p3.numerator,q.p3.denominator))
print('delta  order  | true_direct/rho0 , true_fold/rho0 | formula_direct , formula_fold | ratio_direct , ratio_fold')
for (dl,rc,dt,ft),(dl2,Df,Ff) in zip(res,out):
    dts, fts = sp.simplify(toS(dt)), sp.simplify(toS(ft))
    rd = sp.simplify(dts/Df) if Df!=0 else None
    rf = sp.simplify(fts/Ff) if Ff!=0 else None
    print(dl, rc, '|', dts, ',', fts, '|', Df, ',', Ff, '|', rd, ',', rf)
