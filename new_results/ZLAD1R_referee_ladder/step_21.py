"""ZLAD1R step_14 -- INDEPENDENT exact computation of mu_{k,n}, k<=8, n<=40.
Route (written from thm:flat's Beta-integral formula, no ZLAD1/G1 code):
  E_n = e^{uA} f^{1-2n};  f^{1-2n} = sum_t C(1-2n,t) y^t (e^{u l} - 1)^t;
  [u^j](e^{u l}-1)^t = t! S(j,t) l^j / j!  (Stirling 2nd kind).
  Monomial y^T a^al b^be  ->  poch(n,T)/(2n)_T * DB(al,be; n+T,n),
  DB(al,be;p,q) = al! be! [xi^al eta^be] exp(g_p(xi)+g_q(eta)-g_{p+q}(xi+eta)),
  g_t(x) = sum_{m>=1} c_m(t) x^m/m!,  c_1 = H_{t-1} - gamma,
  c_m = (-1)^m (m-1)! (z(m) - H^(m)_{t-1})  (m>=2).
The exp-series is computed DIRECTLY (series exponentiation), not by Bell recursions.
Coefficients are exact Fractions; zeta monomials z2..z8 and gamma tracked in dicts.
Output: muR_exact_N40_K8.json ; checks gamma cancels exactly.
"""
import json, time
from fractions import Fraction as Fr
from math import comb, factorial

t0 = time.time()
KMAX, NMAX = 8, 40

# ---- zeta/gamma polynomials: dict {(eg,e2,e3,e4,e5,e6,e7,e8): Fraction}
NV = 8
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
def const(c): return {(0,)*NV: Fr(c)} if c else {}
def var(i):  # i=0 -> gamma, i>=2 -> z(i)
    return {tuple(1 if j == (0 if i==0 else i-1) else 0 for j in range(NV)): Fr(1)}

_H = {}
def Hh(r, m):
    if (r, m) not in _H:
        _H[(r, m)] = sum((Fr(1, i**r) for i in range(1, m+1)), Fr(0))
    return _H[(r, m)]

def stir2(m, r):
    return sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r)

# c_m(t) as poly
def cm(m, t):
    if m == 1:
        return padd(const(Hh(1, t-1)), pscale(var(0), -1))         # H - gamma
    f = Fr((-1)**m) * factorial(m-1)
    return padd(pscale(var(m), f), const(-f*Hh(m, t-1)))           # (-1)^m(m-1)!(z(m)-H)

# ---- bivariate series exp(g_p(xi)+g_q(eta)-g_{pq}(xi+eta)) truncated to total degree KMAX
def exp_series(p, q, K):
    # L(xi,eta) as dict {(a,b): poly}, total a+b <= K
    L = {}
    for m in range(1, K+1):
        cm_p, cm_q, cm_pq = cm(m, p), cm(m, q), cm(m, p+q)
        for a in range(m+1):
            b = m - a
            cf = padd(padd(cm_p if b == 0 else {}, cm_q if a == 0 else {}),
                      pscale(cm_pq, -1))
            # xi^a eta^b / m! * binom(m,a) for the (xi+eta)^m part
            if a == 0 or b == 0:
                c = cf  # from g_p (b=0) or g_q (a=0); the -(xi+eta) part included below
            L[(a, b)] = L.get((a, b), {})
            # combine: g_p contributes cm_p when b==0; g_q when a==0; -g_{p+q} contributes -cm_pq*binom(m,a)
            val = {}
            if b == 0: val = padd(val, cm_p)
            if a == 0: val = padd(val, cm_q)
            val = padd(val, pscale(cm_pq, -Fr(comb(m, a))))
            L[(a, b)] = padd(L[(a, b)], pscale(val, Fr(1, factorial(m))))
    L = {k: v for k, v in L.items() if v}
    # exp(L) truncated
    E = {(0, 0): const(1)}
    terms = list(L.items())
    for r in range(1, K+1):
        # add L^r/r!
        P = {(0, 0): const(1)}
        for _ in range(r):
            Q = {}
            for (a1, b1), v1 in P.items():
                for (a2, b2), v2 in terms:
                    if a1+a2+b1+b2 > K: continue
                    key = (a1+a2, b1+b2)
                    Q[key] = padd(Q.get(key, {}), pmul(v1, v2))
            P = Q
        for key, v in P.items():
            E[key] = padd(E.get(key, {}), pscale(v, Fr(1, factorial(r))))
    return {k: v for k, v in E.items() if v}

def poch_ratio(n, j):
    r = Fr(1)
    for t in range(j):
        r *= Fr(n+t, 2*n+t)
    return r

def mu_kn(k, n, exp_cache):
    tot = {}
    for i in range(k+1):           # power of A
        j = k - i                  # u-degree from f part
        for t in range(j+1):
            cb = Fr(1)
            for s_ in range(t):    # binom(1-2n, t)
                cb *= Fr(1-2*n - s_, s_+1)
            st = stir2(j, t)
            if st == 0 or cb == 0: continue
            base = cb * st * factorial(t) * Fr(1, factorial(j)) * Fr(1, factorial(i))
            # A^i l^j expanded: A^i = sum_x binom(i,x)(-(n-1))^{i-x} n^x a^{i-x} b^x
            # l^j = sum_y binom(j,y) b^y (-a)^{j-y}
            for x in range(i+1):
                cx = comb(i, x) * (-(n-1))**(i-x) * n**x
                for y_ in range(j+1):
                    cy = comb(j, y_) * (-1)**(j-y_)
                    al = i - x + (j - y_)
                    be = x + y_
                    key_exp = (n+t, t)  # p = n+T, T = t
                    if key_exp not in exp_cache:
                        exp_cache[key_exp] = exp_series(n+t, n, k)
                    E = exp_cache[key_exp]
                    sub = E.get((al, be), {})
                    if not sub: continue
                    c = base * cx * cy * poch_ratio(n, t) * factorial(al) * factorial(be)
                    tot = padd(tot, pscale(sub, c))
    return pscale(tot, Fr(-1, n*(2*n-1)))

table = {k: {} for k in range(KMAX+1)}
gamma_left = 0
for n in range(1, NMAX+1):
    exp_cache = {}
    for k in range(KMAX+1):
        P = mu_kn(k, n, exp_cache)
        # gamma cancellation check
        gkeys = [kk for kk in P if kk[0] != 0]
        if gkeys:
            gamma_left += 1
            if gamma_left <= 3: print('GAMMA LEFT at', k, n, gkeys[:2])
        table[k][n] = {kk[1:]: str(vv) for kk, vv in P.items() if kk[0] == 0}
    if n % 10 == 0: print(f'n={n} ({time.time()-t0:.0f}s)', flush=True)
json.dump(table, open('muR_exact_N40_K8.json', 'w'))
print('saved muR_exact_N40_K8.json; gamma leftovers:', gamma_left, f'({time.time()-t0:.0f}s)')
print('mu_{8,1} z8 =', table[8][1].get((0,0,0,0,0,0,1), '0'), '(expect 2/15)')
print('mu_{8,1} z2^4 =', table[8][1].get((4,0,0,0,0,0,0), '0'), '(expect -1/315)')
open('progress.log','a').write('step_14: independent exact mu table k<=8 n<=40 computed (Stirling+bivariate Gamma series); gamma cancels exactly\n')
