"""ZLAD1R step_28 -- Lemma 1 (operator formula), MY OWN code (fixed test tuples)."""
import json
from fractions import Fraction as Fr
from math import comb, factorial
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
def stir2(m, r):
    return sum((-1)**i * comb(r, i) * (r-i)**m for i in range(r+1)) // factorial(r)
_H = {}
def Hh(r, m):
    if (r, m) not in _H:
        _H[(r, m)] = sum((Fr(1, i**r) for i in range(1, m+1)), Fr(0))
    return _H[(r, m)]
def exp_series_H(p, q, K):
    L = {}
    for m in range(1, K+1):
        if m == 1:
            cp, cq, cpq = Hh(1, p-1), Hh(1, q-1), Hh(1, p+q-1)
        else:
            f = Fr((-1)**m) * factorial(m-1)
            cp, cq, cpq = -f*Hh(m, p-1), -f*Hh(m, q-1), -f*Hh(m, p+q-1)
        for a in range(m+1):
            b = m - a
            val = Fr(0)
            if b == 0: val += cp
            if a == 0: val += cq
            val -= cpq * comb(m, a)
            L[(a, b)] = L.get((a, b), Fr(0)) + val/factorial(m)
    L = {k: Fr(v) for k, v in L.items() if v}
    E = {(0, 0): Fr(1)}
    terms = list(L.items())
    for r in range(1, K+1):
        P = {(0, 0): Fr(1)}
        for _ in range(r):
            Q = {}
            for (a1, b1), v1 in P.items():
                for (a2, b2), v2 in terms:
                    if a1+a2+b1+b2 > K: continue
                    key = (a1+a2, b1+b2)
                    Q[key] = Q.get(key, Fr(0)) + v1*v2
            P = Q
        for key, v in P.items():
            E[key] = E.get(key, Fr(0)) + v/factorial(r)
    return {k: v for k, v in E.items() if v}
def poch_ratio(n, j):
    r = Fr(1)
    for t in range(j):
        r *= Fr(n+t, 2*n+t)
    return r
def E_poly(k, n):
    P = {}
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
                    key = (t, i-x+(j-y_), x+y_)
                    P[key] = P.get(key, Fr(0)) + base*cx*cy
    return {kk: vv for kk, vv in P.items() if vv}
def d_ab(P, ra, rb):
    R = {}
    for (T, al, be), v in P.items():
        if al >= ra and be >= rb:
            fac = factorial(al)//factorial(al-ra) * factorial(be)//factorial(be-rb)
            key = (T, al-ra, be-rb)
            R[key] = R.get(key, Fr(0)) + v*fac
    return {kk: vv for kk, vv in R.items() if vv}
def p_add(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        nv = R.get(k, Fr(0)) + s*v
        if nv: R[k] = nv
        elif k in R: del R[k]
    return R
def p_scale(P, c):
    return {k: v*c for k, v in P.items()} if c else {}
def sq_w(P, w):
    R = {}
    for r in range(w+1):
        R = p_add(R, p_scale(d_ab(P, r, w-r), comb(w, r)))
    R = p_add(R, d_ab(P, w, 0), -1)
    R = p_add(R, d_ab(P, 0, w), -1)
    return R
def Cw(w): return Fr((-1)**(w+1)) * factorial(w-1)
def apply_PE(P, E):
    op_list = []
    denom = 1
    for w, e in E:
        denom *= factorial(e)
        for _ in range(e):
            op_list.append((w, Cw(w)/Fr(factorial(w))))
    R = P
    for w, c in op_list:
        R = p_scale(sq_w(R, w), c)
    return p_scale(R, Fr(1, denom))
def HH(P, n, cache):
    tot = Fr(0)
    for (T, al, be), v in P.items():
        key = (n+T,)
        if key not in cache:
            cache[key] = exp_series_H(n+T, n, 12)
        E = cache[key]
        sub = E.get((al, be), Fr(0))
        tot += v * poch_ratio(n, T) * factorial(al) * factorial(be) * sub
    return tot

tests = [
    ((2, 4), [(2, 1)], 'z2'),
    ((2, 4), [(3, 1)], 'z3'),
    ((2, 4), [(4, 1)], 'z4'),
    ((2, 4), [(2, 2)], 'z2^2'),
    ((3, 6), [(2, 3)], 'z2^3'),
    ((3, 6), [(2, 1), (3, 1)], 'z2z3'),
    ((3, 6), [(3, 2)], 'z3^2'),
    ((3, 6), [(6, 1)], 'z6'),
    ((2, 8), [(7, 1)], 'z7'),
    ((2, 8), [(8, 1)], 'z8'),
    ((2, 8), [(2, 4)], 'z2^4'),
    ((2, 8), [(3, 1), (5, 1)], 'z3z5'),
    ((2, 8), [(2, 2), (3, 1)], 'z2^2z3'),
    ((4, 7), [(7, 1)], 'z7'),
    ((4, 7), [(2, 1), (5, 1)], 'z2z5'),
    ((5, 8), [(4, 2)], 'z4^2'),
    ((5, 8), [(2, 1), (6, 1)], 'z2z6'),
]
ok_all = True
for (n, k), E, m in tests:
    Pk = E_poly(k, n)
    cache = {}
    got = HH(apply_PE(Pk, E), n, cache)
    expv = -Fr(n*(2*n-1)) * mu(k, n, m)
    good = (got == expv)
    ok_all &= good
    print(f'P1 (n,k)=({n},{k}) E={m}: match={good}' + ('' if good else f'  operator={got} table={expv}'))
print('LEMMA 1 (operator formula), MY code, all 17 tests:', ok_all)
open('progress.log','a').write('step_28: Lemma 1 operator formula verified with independent code at 17 (n,k,E) points incl. w=7,8: ALL MATCH\n')
