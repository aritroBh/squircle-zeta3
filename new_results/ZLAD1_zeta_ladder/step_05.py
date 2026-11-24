"""step_05.py -- verification of the two proof pillars (fixed).
P2: (d_a+d_b)E_n(u) = u E_n(u), exact polynomial check, n=1..3, order u^6.
P1: [z^E]g_{k,n} = [P^E(d_a,d_b) P_k]-Heval, for (n,k)=(2,4) and (3,5), several monomials E.
Heval maps y^j A^i L^l -> poch(n,j)/poch(2n,j) * DH(i,l; n+j) with DH the zeta-FREE directional Bell polys.
"""
import sys, json
sys.path.insert(0, '.')
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
                    res = padd(res, pmul(const(Fr(comb(i-1, r)*comb(l, s))), pmul(kapH(r+1, s, pj), DH(i-1-r, l-s, pj))))
        else:
            res = {}
            for s in range(l):
                res = padd(res, pmul(const(Fr(comb(l-1, s))), pmul(kapH(0, s+1, pj), DH(0, l-1-s, pj))))
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
        tot = padd(tot, pscale(apply_op(P, dA, dL), cv))
    return tot

def apply_op(P, dA, dL):
    R = {}
    for (j, iA, iL), v in P.items():
        if iA >= dA and iL >= dL:
            fac = factorial(iA)//factorial(iA-dA)*factorial(iL)//factorial(iL-dL)
            R[(j, iA-dA, iL-dL)] = R.get((j, iA-dA, iL-dL), 0) + v*fac
    return R

def square(P, w, nv):
    tot = apply_op(P, w, 0)
    tot = padd(tot, pscale(op_on_ab(P, w, 0, nv), -1))
    tot = padd(tot, pscale(op_on_ab(P, 0, w, nv), -1))
    return tot

def Cw(w): return Fr((-1)**(w+1)*factorial(w-1))

def P2_check():
    ok_all = True
    for nv in (1, 2, 3):
        Ek = E_series(nv, 6)
        for k in range(7):
            D = {}
            for (j, iA, iL), v in Ek[k].items():
                if iA > 0:
                    D[(j, iA-1, iL)] = D.get((j, iA-1, iL), 0) + v*iA
            rhs = Ek[k-1] if k >= 1 else {}
            if padd(D, pscale(rhs, -1)):
                ok_all = False
    print('P2: (d_a+d_b)E_n = u E_n for n=1,2,3 to order u^6:', ok_all)
    return ok_all

def P1_check(nv, k, tab=None):
    Ek = E_series(nv, k); Pk = Ek[k]; DH = make_DH(nv)
    if tab is None:
        t = json.load(open('mu_exact_N150_K6.json'))
    else:
        t = tab
    pref = Fr(-nv*(2*nv-1))
    tests = [
        ('z2',   pscale(square(Pk, 2, nv), Cw(2)/Fr(2))),
        ('z2^2', pscale(square(square(Pk, 2, nv), 2, nv), (Cw(2)/Fr(2))**2/2)),
        ('z3',   pscale(square(Pk, 3, nv), Cw(3)/Fr(6))),
        ('z4',   pscale(square(Pk, 4, nv), Cw(4)/Fr(24))),
        ('z2z3', pscale(square(square(Pk, 2, nv), 3, nv), (Cw(2)/Fr(2))*(Cw(3)/Fr(6)))),
        ('z2^3', pscale(square(square(square(Pk, 2, nv), 2, nv), 2, nv), (Cw(2)/Fr(2))**3/6)),
    ]
    ok = True
    for name, Q in tests:
        got = Heval(nv, Q, DH).get((0,)*ZV, Fr(0))
        try:
            expv = t[k][nv].get(name, '0')
        except (KeyError, TypeError):
            expv = t[str(k)][str(nv)].get(name, '0')
        exp = pref*Fr(expv)
        good = (got == exp)
        ok = ok and good
        print(f'  P1 [{name}] g_{{{k},{nv}}}: operator = {got}, table = {exp}, match = {good}')
    print(f'P1 at (n,k)=({nv},{k}):', ok)
    return ok

if __name__ == '__main__':
    ok = P2_check()
    ok = P1_check(2, 4) and ok
    ok = P1_check(3, 5) and ok
    print('ALL PILLARS OK:', ok)
