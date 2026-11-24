
# step_12b.py : PSLQ record with residuals and proper pi->zeta conversion
import mpmath as mp
import sympy as sp
import pickle

c_exact = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}
z2 = sp.zeta(2)
pi2z = {sp.pi**2: 6*z2, sp.pi**4: 36*z2**2, sp.pi**6: 216*z2**3, sp.pi**8: 1296*z2**4}
def to_zeta_basis(e):
    e = sp.expand(e)
    for pw, rep in pi2z.items():
        e = e.subs(pw, rep)
    return sp.expand(e)

def basis_elements(N):
    out = []
    for a in range(N+1):
        for b in range(N//2+1):
            for c in range(N//3+1):
                for d in range(N//5+1):
                    for e in range(N//7+1):
                        if a+2*b+3*c+5*d+7*e <= N:
                            out.append((a,b,c,d,e))
    return out

def monom(m):
    a,b,c,d,e = m
    return sp.log(2)**a * sp.zeta(2)**b * sp.zeta(3)**c * sp.zeta(5)**d * sp.zeta(7)**e

def exact_vector(n, basis):
    e = to_zeta_basis(c_exact[n])
    vec = []
    for m in basis:
        mm = monom(m)
        vec.append(sp.expand(e).coeff(mm))
        e = sp.expand(e - vec[-1]*mm)
    assert sp.expand(e) == 0, f"leftover {e}"
    return vec

for n in [5,6,7,8]:
    basis = basis_elements(n)
    vec = exact_vector(n, basis)
    # primitive integer scaling
    from math import gcd
    from functools import reduce
    dens = [sp.denom(v) for v in vec]
    L = sp.lcm(dens)
    ivec = [int(v*L) for v in vec]
    g = reduce(gcd, [abs(v) for v in ivec if v])
    ivec = [v//g for v in ivec]
    print(f"\nc_{n}: exact integer relation coeffs max={max(abs(v) for v in ivec)} (scale L/g={L//g})")
    for dps in [100, 150]:
        mp.mp.dps = dps
        L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
        target = mp.mpf(sp.N(c_exact[n], dps+10))
        vals = [L2**m[0]*z[2]**m[1]*z[3]**m[2]*z[5]**m[3]*z[7]**m[4] for m in basis]
        rel = mp.pslq([target]+vals, maxcoeff=10**7, maxsteps=10000)
        if rel is None:
            print(f"  dps={dps}: FAILED"); continue
        resid = abs(rel[0]*target + sum(r*v for r,v in zip(rel[1:],vals)))
        # compare rel[1:] with ivec up to scale
        import fractions
        rat = sp.Rational(rel[1:][ivec.index(next(v for v in ivec if v))], next(v for v in ivec if v))
        scaled = [sp.Rational(r)/rat for r in rel[1:]]
        match = all(sp.Rational(r)/rat == v for r,v in zip(rel[1:], ivec))
        print(f"  dps={dps}: residual={mp.nstr(resid,3)}  max|coeff|={max(abs(r) for r in rel)}  matches-exact-vector={match}")
