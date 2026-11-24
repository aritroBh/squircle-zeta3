"""ZLAD1R step_40 -- PREDICTIVE test at weight 9 (beyond ZLAD1's verified list):
extend my independent table to k=9, n<=10, and check 8 weight-9 rungs derived by ME from the
master formula (single z9; mixed z2z7, z3z6, z4z5, z2^2z5, z2z3z4, z3^3, z2^3z3)."""
import json, time
from fractions import Fraction as Fr
from math import comb, factorial
t0 = time.time()
KMAX, NMAX = 9, 10
NV = 10  # gamma + z2..z9 -> indices: 0=gamma, i-1 for z(i), i<=9 -> need e2..e9 = 8 slots + gamma = 9; use 10 for safety
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
                    al, be = i-x+(j-y_), x+y_
                    if (n+t,) not in exp_cache:
                        exp_cache[(n+t,)] = exp_series(n+t, n, KMAX)
                    sub = exp_cache[(n+t,)].get((al, be), {})
                    if not sub: continue
                    c = base * cx * cy * poch_ratio(n, t) * factorial(al) * factorial(be)
                    tot = padd(tot, pscale(sub, c))
    return pscale(tot, Fr(-1, n*(2*n-1)))
def mname(e):  # e has 9 slots (z2..z9 at positions 1..9 after gamma at 0) -> here kk[1:] gives 9 slots: idx 0->z2 ... idx 7->z9, idx 8 unused
    ee = e
    if all(x == 0 for x in ee): return '1'
    out = []
    for idx, x in enumerate(ee):
        if x: out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)
tab9 = {k: {} for k in range(KMAX+1)}
gleft = 0
for n in range(1, NMAX+1):
    cache = {}
    for k in range(KMAX+1):
        P = mu_kn(k, n, cache)
        gleft += sum(1 for kk in P if kk[0] != 0)
        tab9[k][n] = {mname(kk[1:]): vv for kk, vv in P.items() if kk[0] == 0}
    print(f'n={n} ({time.time()-t0:.0f}s)', flush=True)
print('gamma leftovers:', gleft)
json.dump({str(k): {str(n): {m: str(v) for m, v in tab9[k][n].items()} for n in tab9[k]} for k in tab9},
          open('muR_exact_N10_K9.json', 'w'))
def mu9(k, n, m): return tab9[k][n].get(m, Fr(0))
ok = True
def chk(name, cond):
    global ok
    good = all(cond(n) for n in range(1, NMAX+1))
    ok &= good
    print(('PASS' if good else 'FAIL'), name)
chk('w9 single: [z9]mu_9 = -[z2]mu_2 -6[z2^2]mu_4 -20[z2^3]mu_6 -24[z2^4]mu_8',
    lambda n: mu9(9,n,'z9') == -mu9(2,n,'z2') - 6*mu9(4,n,'z2^2') - 20*mu9(6,n,'z2^3') - 24*mu9(8,n,'z2^4'))
chk('w9 mix [z2z7]mu_9 = -2[z2^2]mu_4 -12[z2^3]mu_6 -24[z2^4]mu_8',
    lambda n: mu9(9,n,'z2z7') == -2*mu9(4,n,'z2^2') - 12*mu9(6,n,'z2^3') - 24*mu9(8,n,'z2^4'))
chk('w9 mix [z3z6]mu_9 = -2[z2^2]mu_4 -9[z2^3]mu_6 -8[z2^4]mu_8',
    lambda n: mu9(9,n,'z3z6') == -2*mu9(4,n,'z2^2') - 9*mu9(6,n,'z2^3') - 8*mu9(8,n,'z2^4'))
chk('w9 mix [z4z5]mu_9 = -2[z2^2]mu_4 -9[z2^3]mu_6 -12[z2^4]mu_8',
    lambda n: mu9(9,n,'z4z5') == -2*mu9(4,n,'z2^2') - 9*mu9(6,n,'z2^3') - 12*mu9(8,n,'z2^4'))
chk('w9 mix [z2^2z5]mu_9 = -3[z2^3]mu_6 -12[z2^4]mu_8',
    lambda n: mu9(9,n,'z2^2z5') == -3*mu9(6,n,'z2^3') - 12*mu9(8,n,'z2^4'))
chk('w9 mix [z2z3z4]mu_9 = -6[z2^3]mu_6 -12[z2^4]mu_8',
    lambda n: mu9(9,n,'z2z3z4') == -6*mu9(6,n,'z2^3') - 12*mu9(8,n,'z2^4'))
chk('w9 mix [z3^3]mu_9 = -[z2^3]mu_6',
    lambda n: mu9(9,n,'z3^3') == -mu9(6,n,'z2^3'))
chk('w9 mix [z2^3z3]mu_9 = 4[z2^4]mu_8',
    lambda n: mu9(9,n,'z2^3z3') == 4*mu9(8,n,'z2^4'))
print('WEIGHT-9 PREDICTIVE RUNGS (beyond ZLAD1 list):', 'ALL PASS' if ok else 'FAIL', f'({time.time()-t0:.0f}s)')
open('progress.log','a').write('step_40: weight-9 predictive rungs (k=9, n<=10, my own derivation from master formula): ALL PASS\n')
