import sys, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
from math import comb, factorial
from step_02 import E_series, poch_ratio, const, padd, pmul, pscale, H, ZV

def psiH(m, t):
    if m == 1: return const(H(1, t-1))
    f = Fr((-1)**m)
    for i in range(1, m): f *= i
    return const(-f*H(m, t-1))
def make_DH(nv):
    def kapH(r, s, pj):
        p = nv+pj; q = nv; m = r+s
        out = pmul(const(Fr((-(nv-1))**r * (-1)**s)), psiH(m, p))
        out = padd(out, pmul(const(Fr(nv**r)), psiH(m, q)))
        if s == 0: out = padd(out, pscale(psiH(m, p+q), -1))
        return out
    cache = {}
    def DH(i, l, pj):
        key = (i, l, pj)
        if key in cache: return cache[key]
        if i == 0 and l == 0: res = const(1)
        elif i > 0:
            res = {}
            for r in range(i):
                for s in range(l+1):
                    res = padd(res, pmul(const(Fr(comb(i-1,r)*comb(l,s))), pmul(kapH(r+1,s,pj), DH(i-1-r,l-s,pj))))
        else:
            res = {}
            for s in range(l):
                res = padd(res, pmul(const(Fr(comb(l-1,s))), pmul(kapH(0,s+1,pj), DH(0,l-1-s,pj))))
        cache[key] = res
        return res
    return DH
def Heval(nv, P, DH):
    tot = {}
    for (j, iA, iL), cv in P.items():
        tot = padd(tot, pmul(const(cv*poch_ratio(nv, j)), DH(iA, iL, j)))
    return tot
def op_on_ab(P, ra, sb, nv):
    ops = {}
    for x in range(ra+1):
        for y2 in range(sb+1):
            cA = comb(ra, x) * (-(nv-1))**(ra-x) * (-1)**x
            cB = comb(sb, y2) * nv**(sb-y2)
            key = (ra-x+sb-y2, x+y2)
            ops[key] = ops.get(key, 0) + cA*cB
    tot = {}
    for (dA, dL), cv in ops.items():
        R = {}
        for (j, iA, iL), v in P.items():
            if iA >= dA and iL >= dL:
                key = (j, iA-dA, iL-dL)
                fac = factorial(iA)//factorial(iA-dA)*factorial(iL)//factorial(iL-dL)
                R[key] = R.get(key, 0) + v*fac
        tot = padd(tot, pscale(R, cv))
    return tot
def square(P, w, nv):
    tot = {}
    # dA^w
    R = {}
    for (j,iA,iL),v in P.items():
        if iA >= w:
            fac = factorial(iA)//factorial(iA-w)
            R[(j,iA-w,iL)] = R.get((j,iA-w,iL),0)+v*fac
    tot = padd(tot, R)
    tot = padd(tot, pscale(op_on_ab(P, w, 0, nv), -1))
    tot = padd(tot, pscale(op_on_ab(P, 0, w, nv), -1))
    return tot

nv, k = 2, 4
Ek = E_series(nv, k); P4 = Ek[k]; DH = make_DH(nv)
t = json.load(open('mu_exact_N150_K6.json'))
pref = Fr(-nv*(2*nv-1))
def Cw(w): return Fr((-1)**(w+1)*factorial(w-1))
# rational part:
print('rational part g_{4,2}:', Heval(nv, P4, DH)[(0,)*ZV], 'exp:', pref*Fr(t['4']['2']['1']))
tests = [
 ('z2',   pscale(square(P4, 2, nv), Cw(2)/Fr(2))),
 ('z2^2', pscale(square(square(P4, 2, nv), 2, nv), (Cw(2)/Fr(2))**2 / 2)),
 ('z3',   pscale(square(P4, 3, nv), Cw(3)/Fr(6))),
 ('z4',   pscale(square(P4, 4, nv), Cw(4)/Fr(24))),
 ('z2z3', pscale(square(square(P4, 2, nv), 3, nv), (Cw(2)/Fr(2))*(Cw(3)/Fr(6)))),
]
for name, Q in tests:
    got = Heval(nv, Q, DH).get((0,)*ZV, Fr(0))
    expv = t['4']['2'].get(name)
    exp = pref*Fr(expv) if expv is not None else None
    print(f'[{name}]g_{{4,2}}: operator = {got}  expected = {exp}  match: {got == exp}')
