
# step_12c.py : PSLQ record, fixed
import mpmath as mp
import sympy as sp
import pickle
from math import gcd
from functools import reduce

c_exact = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}
Z2 = sp.Symbol('Z2'); dpi = sp.Symbol('dpi')
def to_zeta_basis(e):
    return sp.expand(sp.expand(e).subs(sp.pi, dpi).subs(dpi**2, 6*Z2))

def basis_elements(N):
    out = []
    for a in range(N+1):
        for b in range(N//2+1):
            for c in range(N//3+1):
                for dd in range(N//5+1):
                    for e in range(N//7+1):
                        if a+2*b+3*c+5*dd+7*e <= N: out.append((a,b,c,dd,e))
    return out
def monom(m):
    a,b,c,dd,e = m
    return sp.log(2)**a * Z2**b * sp.zeta(3)**c * sp.zeta(5)**dd * sp.zeta(7)**e

for n in [5,6,7,8]:
    basis = basis_elements(n)
    e = to_zeta_basis(c_exact[n])
    vec = []
    for m in basis:
        mm = monom(m)
        ex = sp.expand(e)
        cf = ex.coeff(Z2,m[1]).coeff(sp.log(2),m[0]).coeff(sp.zeta(3),m[2]).coeff(sp.zeta(5),m[3]).coeff(sp.zeta(7),m[4])
        vec.append(cf); e = sp.expand(e - cf*mm)
    assert sp.expand(e)==0, f"leftover {e}"
    dens = [sp.denom(v) for v in vec]; L = sp.lcm(dens)
    ivec = [int(v*L) for v in vec]; g = reduce(gcd,[abs(v) for v in ivec if v]); ivec=[v//g for v in ivec]
    print(f"\nc_{n}: exact primitive integer vector max|coeff| = {max(abs(v) for v in ivec)}")
    for dps in [100,150]:
        mp.mp.dps = dps
        L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
        target = mp.mpf(sp.N(c_exact[n], dps+10))
        vals = [L2**m[0]*z[2]**m[1]*z[3]**m[2]*z[5]**m[3]*z[7]**m[4] for m in basis]
        rel = mp.pslq([target]+vals, maxcoeff=10**7, maxsteps=10000)
        resid = abs(rel[0]*target + sum(r*v for r,v in zip(rel[1:],vals)))
        i0 = next(i for i,v in enumerate(ivec) if v)
        rat = sp.Rational(rel[1:][i0], ivec[i0])
        match = all(sp.Rational(r) == rat*v for r,v in zip(rel[1:],ivec))
        print(f"  dps={dps}: resid={mp.nstr(resid,3)} max|c|={max(abs(r) for r in rel)} match={match}")

print("\n-- incomplete-basis controls (dps=100) --")
for n, drop, tag in [(5,(0,0,1,0,0),'no zeta(3)'), (6,(0,0,0,1,0),'no zeta(5)'), (8,(0,0,0,0,1),'no zeta(7)')]:
    mp.mp.dps = 100
    L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
    target = mp.mpf(sp.N(c_exact[n], 110))
    basis = basis_elements(n)
    basis = [m for m in basis if not (m[2]>0 and drop[2]) and not (m[3]>0 and drop[3]) and not (m[4]>0 and drop[4])]
    vals = [L2**m[0]*z[2]**m[1]*z[3]**m[2]*z[5]**m[3]*z[7]**m[4] for m in basis]
    rel = mp.pslq([target]+vals, maxcoeff=10**6, maxsteps=10000)
    print(f"c_{n} {tag}: {'FAILED (correctly - element is required)' if rel is None else 'found maxc='+str(max(abs(r) for r in rel))+' resid='+mp.nstr(abs(rel[0]*target+sum(rr*v for rr,v in zip(rel[1:],vals))),3)}")
