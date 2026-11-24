
lib = '''"""Shared helpers for CHN2 Part C."""
from fractions import Fraction as F
import mzv_exact as M
import engine5

def peval5(p, a):
    return sum(c*a[0]**e[0]*a[1]**e[1]*a[2]**e[2]*a[3]**e[3]*a[4]**e[4] for e,c in p.items())

def gamprod_series_point(coeffs, Wmax):
    q=[{}]*(Wmax+1)
    for m in range(2,Wmax+1):
        s=sum(F(e)*ell**m for ell,e in coeffs)
        q[m]={(m,): F((-1)**m,m)*s}
    E=[{} for _ in range(Wmax+1)]; E[0]={():F(1)}
    for W in range(1,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for z,c in q[i].items():
                for z2,c2 in E[W-i].items():
                    z3=tuple(sorted(z+z2))
                    acc[z3]=acc.get(z3,F(0))+F(i,W)*c*c2
        E[W]=acc
    return E

def m0_coords(a, Wmax):
    """m=0 (Gamma-product prefactor) contribution to c_W coords at rational point a (tuple of 5 rationals)."""
    a1,a2,a3,a4,a5 = [F(x) for x in a]
    S4 = a1+a2+a3+a4; S4r = a5+a4+a3+a2
    s1 = gamprod_series_point([(a2+a4,1),(S4,-1),(a1+a3,1)], Wmax)
    s2 = gamprod_series_point([(a2+a4,1),(S4r,-1),(a5+a3,1)], Wmax)
    f1 = a2/(a2+a4); f2 = a4/(a2+a4)
    out={}
    for W in range(2,Wmax+1):
        T=engine5.table(W); d=T.nullity
        acc=[F(0)]*d
        for s_,f_ in ((s1,f1),(s2,f2)):
            for zmon,c in s_[W].items():
                c=c*f_
                wc=M.wordcombo(tuple((k,) for k in zmon))
                for w,cc in wc.items():
                    rw=T.reduce_word(w)
                    for b in range(d): acc[b]+=c*cc*rw[b]
        out[W]=acc
    return out

def full_coords_at(a, Wmax, forms=None):
    """engine regular (from forms) + m0 at point a -> {W: coord list}."""
    if forms is None:
        forms = engine5.reduced_forms(Wmax)[0]
    out={}
    for W in range(2,Wmax+1):
        T=engine5.table(W); d=T.nullity
        if W>=3:
            reg=[peval5(forms[W][b], a) for b in range(d)]
        else:
            reg=[F(0)]*d
        m0=m0_coords(a, W)[W]
        out[W]=[reg[b]+m0[b] for b in range(d)]
    return out
'''
open('kimilib.py','w').write(lib)
import importlib
import kimilib, engine5
importlib.reload(engine5); importlib.reload(kimilib)
from fractions import Fraction as F
truth = {3:[F(10)], 4:[F(-19)], 5:[F(102),F(-36)], 6:[F(-159),F(68)]}
m0 = kimilib.m0_coords((0,1,1,1,1), 6)
forms6, eng6 = engine5.reduced_forms(6)
for W in range(3,7):
    reg = [kimilib.peval5(forms6[W][b], (0,1,1,1,1)) for b in range(len(truth[W]))]
    m0v = m0[W]
    treg = [truth[W][b]-m0v[b] for b in range(len(truth[W]))]
    print(f"W={W}: truth={truth[W]} m0={m0v} => true regular={treg}  my regular={reg}")
for W in [4,5]:
    print(f"--- Areg[{W}] ---")
    for w,p in eng6['Areg'][W].items():
        print("  ", w, p)
