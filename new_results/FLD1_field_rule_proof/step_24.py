WD='<project>/research/08_theory/FLD1_field_rule_proof'
code = r'''"""step22_musymbolic.py -- SYMBOLIC-n closed formula for mu_{k,n}, general k, and its Z-sum canonical form.

mu_{k,n} = -(1/(n(2n-1))) sum_{(j,ia,lb)} e^{(k)}_{j,ia,lb}(n) (n)_j/(2n)_j D_{ia,lb}(n+j,n),
[u^k]E_n = sum e y^j a^ia b^lb via the Stirling formula (binom(1-2n,m) is a polynomial in n; fast).
D = moment polynomials of the cumulants; psi at n+j, n, 2n+j reduces to
   h_r = H^(r)_{n-1}, g_r = H^(r)_{2n-1}, z_r = zeta(r),  with Q(n) coefficients.
Then to_zsum: h_r -> 2^{r-1}(Z_r(2n)+Z_{-r}(2n)) - 1/n^r,  g_r -> Z_r(2n) - 1/(2n)^r,
products -> stuffle.  Output: Zsum[k][zeta_mono] = {word: coef in Q(n)}.
Verified against the fast numeric table (step_21) for n<=60.
"""
import sys, json, pickle
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp
from sympy import QQ
from sympy.polys.fields import field
from sympy.polys.rings import ring
sys.path.insert(0, '.')
import zsum
from zsum import K, nn, ONE, ZERO

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
n = sp.Symbol('n')

def S2(j, m):
    if m > j or m < 0: return 0
    if j == 0: return 1 if m == 0 else 0
    return sum(((-1)**(m-t))*comb(m, t)*t**j for t in range(m+1))//factorial(m)

# ---- [u^k]E_n as dict {(j,ia,lb): sympy poly in n}
def Ek_coeffs(k):
    out = {}
    def add(j,ia,lb,c):
        out[(j,ia,lb)] = sp.expand(out.get((j,ia,lb),sp.Integer(0)) + c)
    for i in range(k+1):
        j = k-i
        fi = Fr(1, factorial(i)); fj = Fr(1, factorial(j))
        # A^i = (-(n-1)a + n b)^i
        for ta in range(i+1):
            tb = i-ta
            ca = Fr(comb(i,ta))
            coefA = ca * (-(n-1))**ta * n**tb   # sympy poly in n
            for m in range(j+1):
                bc = Fr(1)
                for t in range(m): bc *= Fr(1, t+1)
                bcm = sp.expand(sp.Rational(1,1)*sp.prod([sp.Rational(1,t+1)*(1-2*n-t) for t in range(m)])) if m>0 else sp.Integer(1)
                cm = bcm * factorial(m) * S2(j, m)
                for tal in range(j+1):
                    cl = Fr(comb(j,tal))*Fr(-1)**tal
                    # a-exp = ta+tal, b-exp = tb+(j-tal)
                    add(m, ta+tal, tb+(j-tal), sp.Rational(fi.numerator, fi.denominator)*sp.Rational(fj.numerator,fj.denominator)*coefA*cm*sp.Rational(cl.numerator,cl.denominator))
    return out

# ---- ring in h1..hK, g1..gK, z2..zK over K=Q(n)
Hv = [f'h{r}' for r in range(1, KMAX+1)]
Gv = [f'g{r}' for r in range(1, KMAX+1)]
Zv = [f'z{r}' for r in range(2, KMAX+1)]
R = ring(', '.join(Hv+Gv+Zv), K)[0]
gens = R.gens
h = {r: gens[r-1] for r in range(1, KMAX+1)}
g = {r: gens[KMAX + r-1] for r in range(1, KMAX+1)}
z = {r: gens[2*KMAX + r-2] for r in range(2, KMAX+1)}
def toK(expr):
    expr = sp.expand(expr)
    P = sp.Poly(expr, n)
    out = K(0)
    for (e,), c in P.terms():
        c = sp.Rational(c); out += K(QQ(int(c.p), int(c.q)))*nn**e
    return out

def shift_h(r, j):   # H^(r)_{n+j-1} = h_r + sum_{t<j} 1/(n+t)^r
    return h[r] + sum((K(1)/(nn+t)**r for t in range(j)), K(0))
def shift_g(r, j):   # H^(r)_{2n+j-1}
    return g[r] + sum((K(1)/(2*nn+t)**r for t in range(j)), K(0))

def psi(m, which, j):
    base = nn if which=='n' else 2*nn
    H = h if which=='n' else g
    sh = shift_h if which=='n' else shift_g
    if m == 0:
        return sh(1, j)
    c = (-1)**(m+1)*factorial(m)
    # (-1)^{m+1} m! (zeta(m+1) - H^(m+1)_{.-1})
    harm = sh(m+1, j) if (m+1) <= KMAX else None
    out = K(0)
    if (m+1) in z: out += c*z[m+1]
    if harm is not None: out += -c*harm
    return out

def cumulant(r, m, j):
    # kap_{r,m} at p=n+j, q=n (so p+q = 2n+j)
    out = -psi(r+m-1, 'g', j)
    if m == 0: out = out + psi(r+m-1, 'h', j)   # psi(p) = psi(n+j)
    if r == 0: out = out + psi(r+m-1, 'h', 0)   # psi(q) = psi(n)
    return out

def moments(j):
    D = {(0,0): R(1)}
    for tot in range(1, KMAX+1):
        for I in range(tot, -1, -1):
            L = tot-I
            v = R(0)
            if I>0:
                for r in range(1, I+1):
                    for m in range(0, L+1):
                        v += comb(I-1,r-1)*comb(L,m)*cumulant(r,m,j)*D[(I-r,L-m)]
            else:
                for m in range(1, L+1):
                    v += comb(L-1,m-1)*cumulant(0,m,j)*D[(0,L-m)]
            D[(I,L)] = v
    return D

def mu_formula(k):
    tot = R(0)
    Dc = {}
    for (j,ia,lb), co in Ek_coeffs(k).items():
        if j not in Dc: Dc[j] = moments(j)
        pr = K(1)
        for t in range(j): pr *= (nn+t)/(2*nn+t)
        tot += toK(co)*pr*Dc[j][(ia,lb)]
    return tot*(-K(1)/(nn*(2*nn-1)))

# ---- numeric eval of the ring formula at n=nv, compared to step21 table
import json as _json
def zeta_of_exp(e):
    # e is exponent tuple over (h's, g's, z's); extract zeta monomial
    return tuple(e[2*KMAX:])

if __name__ == '__main__':
    table = _json.load(open('mu_exact_K5_N100.json'))
    names = Hv+Gv+Zv
    allok = True
    Zstore = {}
    for k in range(KMAX+1):
        F = mu_formula(k)
        # verify against fast table
        def evalF(nv):
            import mpmath  # not needed; use fractions
        # build evaluator: substitute h_r=H(r,nv-1), g_r=H(r,2nv-1), z_r=zeta (rational? no).
        # Instead compare RATIONAL part and each zeta-monomial separately with step21 (which is exact in Q[z]).
        # Evaluate coefficients of each zeta-monomial as rational numbers by setting h,g to harmonic values, z symbolically tracked.
        def Hf(r,m): return sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
        def evalK(c,nv): return Fr(int(c.numer(nv)),1)/Fr(int(c.denom(nv)),1)
        okk=True
        for nv in range(1,61):
            # for each monomial, coefficient; substitute h,g numeric, keep z-exp as key
            acc={}
            for mon, cf in F.terms():
                cval = evalK(cf, nv)
                hv=1
                for r in range(1,KMAX+1):
                    hv *= Hf(r,nv-1)**mon[r-1]
                    hv *= Hf(r,2*nv-1)**mon[KMAX+r-1]
                ze = tuple(mon[2*KMAX:])
                # pad ze to length 4 (z2..z5)
                ze = tuple(list(ze)+[0]*(4-len(ze)))[:4]
                acc[ze]=acc.get(ze,Fr(0))+cval*hv
            ref = table[str(k)][str(nv)]
            for kk in set(list(acc)+list(ref)):
                if Fr(acc.get(kk,'0'))!=Fr(ref.get(kk,'0')):
                    okk=False; allok=False; print('MISMATCH',k,nv,kk,acc.get(kk),ref.get(kk))
        print(f'k={k}: symbolic formula == fast numeric table (n<=60): {okk};  #monomials={len(F.terms())}')
        sys.stdout.flush()
        Zstore[k]=F
    print('ALL k symbolic==numeric:', allok)
    pickle.dump({k:[(mon, str(cf.as_expr())) for mon,cf in Zstore[k].terms()] for k in Zstore}, open('mu_symbolic_K%d.pkl'%KMAX,'wb'))
    print('saved mu_symbolic_K%d.pkl'%KMAX)
'''
open(WD+'/step_22_musymbolic.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_22_musymbolic.py','5'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-3000:] if r.returncode else '')