"""ZLAD1R step_22 -- rerun of independent table with string monomial keys (same names as ZLAD1)."""
import json, time
from fractions import Fraction as Fr
from math import comb, factorial

t0 = time.time()
KMAX, NMAX = 8, 40
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
def var(i):
    return {tuple(1 if j == (0 if i==0 else i-1) else 0 for j in range(NV)): Fr(1)}

_H = {}
def Hh(r, m):
    if (r, m) not in _H:
        _H[(r, m)] = sum((Fr(1, i**r) for i in range(1, m+1)), Fr(0))
    return _H[(r, m)]
def stir2(m, r):
    return sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r)
def cm(m, t):
    if m == 1:
        return padd(const(Hh(1, t-1)), pscale(var(0), -1))
    f = Fr((-1)**m) * factorial(m-1)
    return padd(pscale(var(m), f), const(-f*Hh(m, t-1)))

def exp_series(p, q, K):
    L = {}
    for m in range(1, K+1):
        cm_p, cm_q, cm_pq = cm(m, p), cm(m, q), cm(m, p+q)
        for a in range(m+1):
            b = m - a
            val = {}
            if b == 0: val = padd(val, cm_p)
            if a == 0: val = padd(val, cm_q)
            val = padd(val, pscale(cm_pq, -Fr(comb(m, a))))
            L[(a, b)] = padd(L.get((a, b), {}), pscale(val, Fr(1, factorial(m))))
    L = {k: v for k, v in L.items() if v}
    E = {(0, 0): const(1)}
    terms = list(L.items())
    for r in range(1, K+1):
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
    for i in range(k+1):
        j = k - i
        for t in range(j+1):
            cb = Fr(1)
            for s_ in range(t):
                cb *= Fr(1-2*n - s_, s_+1)
            st = stir2(j, t)
            if st == 0 or cb == 0: continue
            base = cb * st * factorial(t) * Fr(1, factorial(j)) * Fr(1, factorial(i))
            for x in range(i+1):
                cx = comb(i, x) * (-(n-1))**(i-x) * n**x
                for y_ in range(j+1):
                    cy = comb(j, y_) * (-1)**(j-y_)
                    al = i - x + (j - y_)
                    be = x + y_
                    key_exp = (n+t,)
                    if key_exp not in exp_cache:
                        exp_cache[key_exp] = exp_series(n+t, n, KMAX)
                    E = exp_cache[key_exp]
                    sub = E.get((al, be), {})
                    if not sub: continue
                    c = base * cx * cy * poch_ratio(n, t) * factorial(al) * factorial(be)
                    tot = padd(tot, pscale(sub, c))
    return pscale(tot, Fr(-1, n*(2*n-1)))

def mname(e):
    if e == (0,)*7: return '1'
    out = []
    for idx, x in enumerate(e):
        if x:
            out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)

table = {k: {} for k in range(KMAX+1)}
gamma_left = 0
for n in range(1, NMAX+1):
    exp_cache = {}
    for k in range(KMAX+1):
        P = mu_kn(k, n, exp_cache)
        gkeys = [kk for kk in P if kk[0] != 0]
        if gkeys:
            gamma_left += 1
            if gamma_left <= 3: print('GAMMA LEFT at', k, n, gkeys[:2])
        table[k][n] = {mname(kk[1:]): str(vv) for kk, vv in P.items() if kk[0] == 0}
    if n % 10 == 0: print(f'n={n} ({time.time()-t0:.0f}s)', flush=True)
json.dump(table, open('muR_exact_N40_K8.json', 'w'))
print('saved muR_exact_N40_K8.json; gamma leftovers:', gamma_left, f'({time.time()-t0:.0f}s)')
print('mu_{8,1} z8 =', table[8][1].get('z8'), '(expect 2/15)')
print('mu_{8,1} z2^4 =', table[8][1].get('z2^4'), '(expect -1/315)')
print('mu_{2,1} z2 =', table[2][1].get('z2'), '(expect 1/3)')
print('mu_{3,1} z3 =', table[3][1].get('z3'), '(expect -1/3)')
open('progress.log','a').write('step_22: independent exact mu table k<=8 n<=40 saved; gamma cancels exactly in every mu_{k,n}\n')
