code = r'''"""step_02.py -- EXACT mu_{k,n} for k<=6, n<=NMAX, split by zeta monomial.

Method: per integer n, build [u^k]E_n as a polynomial in (y,A,L) with Fraction coeffs
(E_n = exp(uA) * exp((1-2n) * sum_{m>=1} (uL)^m/m! * L_m(y)), L_m(y) = sum_r (-1)^{r+1}(r-1)! S(m,r) y^r).
Then integrate: y^j A^i L^l -> poch(n,j)/poch(2n,j) * DB(i,l; n+j),
DB = mixed Bell polynomials in directional cumulants
  kap_{r,s}(p,q) = (-(n-1))^r (-1)^s psi^{(r+s-1)}(p) + n^r psi^{(r+s-1)}(q) - [s==0] psi^{(r+s-1)}(p+q),
psi^{(m-1)}(t): m=1 -> H^{(1)}_{t-1}; m>=2 -> (-1)^m (m-1)! (z(m) - H^{(m)}_{t-1}).
Everything exact in Q[z2,z3,z4,z5,z6].
Output: mu_exact_N{N}_K6.json  {k: {n: {monom: coeff}}} with monom keys like '1','z2','z3','z2^2','z2z3','z3^2',...
"""
import sys, json, time
from fractions import Fraction as Fr
from math import comb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 150
KMAX = 6
t0 = time.time()

# ---------- polynomial in zeta vars: dict {(e2,e3,e4,e5,e6): Fraction} ----------
ZV = 5  # number of zeta vars z2..z6
def padd(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        nv = R.get(k, 0) + s*v
        if nv: R[k] = nv
        elif k in R: del R[k]
    return R
def pmul(P, Q):
    R = {}
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            k = tuple(x+y for x, y in zip(k1, k2))
            R[k] = R.get(k, 0) + v1*v2
    return {k: v for k, v in R.items() if v}
def pscale(P, c):
    return {k: v*c for k, v in P.items()} if c else {}
def const(c): return {(0,)*ZV: Fr(c)} if c else {}
def zeta_var(w):  # z(w) as polynomial
    return {tuple(1 if i == w-2 else 0 for i in range(ZV)): Fr(1)}

# harmonic numbers (cached)
_H = {}
def H(r, m):
    key = (r, m)
    if key not in _H:
        _H[key] = sum((Fr(1, i**r) for i in range(1, m+1)), Fr(0))
    return _H[key]

def psi(m, t):
    """psi^{(m-1)}(t) for integer t>=1, m>=1, in Q[z2..z6] (gamma-free by cancellation)."""
    if m == 1:
        return const(H(1, t-1))
    sgn = -1 if m % 2 else 1   # (-1)^m
    return padd(pscale(zeta_var(m), sgn*Fr(comb(m-1,1)*1 if False else 1)*1), pscale(const(H(m, t-1)), -sgn*1))

# fix psi: (-1)^m (m-1)! (z(m) - H)
def psi(m, t):
    if m == 1:
        return const(H(1, t-1))
    f = Fr((-1)**m)
    for i in range(1, m): f *= i  # (-1)^m (m-1)!
    return padd(pscale(zeta_var(m), f), const(-f*H(m, t-1)))

# ---------- E-series per n : poly in (y,A,L) dict {(j,iA,iL): Fraction} per u-degree ----------
def stirling2(m, r):
    return comb(m, r) * sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r) if False else _st2(m, r)
from math import factorial
def _st2(m, r):
    return sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r)

def poly_mul(P, Q):  # multiply polys in (y,A,L)
    R = {}
    for (j1, a1, l1), v1 in P.items():
        for (j2, a2, l2), v2 in Q.items():
            k = (j1+j2, a1+a2, l1+l2)
            R[k] = R.get(k, 0) + v1*v2
    return {k: v for k, v in R.items() if v}

def E_series(nv, K):
    # B(u) = (1-2n) sum_{m=1}^K (uL)^m/m! L_m(y)  as list [B_1..B_K] of polys (coeff of u^m)
    B = [None]*(K+1)
    for m in range(1, K+1):
        P = {}
        for r in range(1, m+1):
            c = (-1)**(r+1) * factorial(r-1) * _st2(m, r)
            if c == 0: continue
            # y^r L^m with coeff (1-2n)*c/m!
            key = (r, 0, m)
            P[key] = P.get(key, 0) + (1-2*nv)*Fr(c, factorial(m))
        B[m] = P
    # exp(B) truncated to degree K
    E = [None]*(K+1)
    E[0] = {(0, 0, 0): Fr(1)}
    for k in range(1, K+1):
        P = {}
        for m in range(1, k+1):
            T = poly_mul(B[m], E[k-m])
            for key, v in T.items():
                P[key] = P.get(key, 0) + v*m
        E[k] = {key: Fr(v, k) for key, v in P.items()}  # E_k = (1/k) sum m B_m E_{k-m}
    # multiply by exp(uA): add A-power
    F = [None]*(K+1)
    for k in range(K+1):
        P = {}
        for i in range(k+1):
            src = E[k-i]
            for (j, a, l), v in src.items():
                key = (j, a+i, l)
                P[key] = P.get(key, 0) + v/factorial(i)
        F[k] = P
    return F

# ---------- directional cumulants and Bell DB(i,l) at (p, q=n) ----------
def make_DB(nv, K):
    # kap[(r,s,pj)] for p = n+j, q = n ; return function DB(i,l,pj)
    def kap(r, s, pj):
        p = nv + pj; q = nv; m = r+s
        out = pmul(const(Fr((-(nv-1))**r * (-1)**s)), psi(m, p))
        out = padd(out, pmul(const(Fr(nv**r)), psi(m, q)))
        if s == 0:
            out = padd(out, pscale(psi(m, p+q), -1))
        return out
    cache = {}
    def DB(i, l, pj):
        key = (i, l, pj)
        if key in cache: return cache[key]
        if i == 0 and l == 0:
            res = const(1)
        elif i > 0:
            res = {}
            for r in range(i):
                for s in range(l+1):
                    T = pmul(const(Fr(comb(i-1, r)*comb(l, s))), pmul(kap(r+1, s, pj), DB(i-1-r, l-s, pj)))
                    res = padd(res, T)
        else:
            res = {}
            for s in range(l):
                T = pmul(const(Fr(comb(l-1, s))), pmul(kap(0, s+1, pj), DB(0, l-1-s, pj)))
                res = padd(res, T)
        cache[key] = res
        return res
    return DB

def poch_ratio(nv, j):
    r = Fr(1)
    for t in range(j):
        r *= Fr(nv+t, 2*nv+t)
    return r

MONO_NAMES = {}
def mono_name(e):
    if e == (0,)*ZV: return '1'
    parts = []
    for i, x in enumerate(e):
        if x:
            parts.append(f'z{i+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(parts)

def mu_all(nv, K, Ek, DB):
    out = {}
    for k in range(K+1):
        tot = {}
        for (j, i, l), cv in Ek[k].items():
            T = pmul(const(cv*poch_ratio(nv, j)), DB(i, l, j))
            tot = padd(tot, T)
        tot = pscale(tot, Fr(-1, nv*(2*nv-1)))
        out[k] = {mono_name(e): str(v) for e, v in tot.items()}
    return out

if __name__ == '__main__':
    table = {k: {} for k in range(KMAX+1)}
    for nv in range(1, NMAX+1):
        Ek = E_series(nv, KMAX)
        DB = make_DB(nv, KMAX)
        row = mu_all(nv, KMAX, Ek, DB)
        for k in range(KMAX+1):
            table[k][nv] = row[k]
        if nv % 25 == 0:
            print(f'n={nv} done ({time.time()-t0:.0f}s)', flush=True)
    fn = f'mu_exact_N{NMAX}_K6.json'
    json.dump(table, open(fn, 'w'))
    print('saved', fn, f'({time.time()-t0:.0f}s)')
    # quick canonical checks
    ok2 = all(Fr(table[2][str(n)].get('z2', '0'))*(4*n*n-1) == 1 for n in range(1, NMAX+1))
    ok3 = all(Fr(table[3][str(n)].get('z3', '0'))*(4*n*n-1) == -1 for n in range(1, NMAX+1))
    print('canonical [z2]mu_2=1/(4n^2-1):', ok2, ' [z3]mu_3=-1/(4n^2-1):', ok3)
    for n in [1, 2, 3]:
        print('mu_5,%d:' % n, table[5][str(n)])
        print('mu_6,%d:' % n, table[6][str(n)])
'''
open('step_02.py','w').write(code)
print('written', len(code))
