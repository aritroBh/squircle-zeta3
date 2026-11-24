code = r'''"""step_05.py -- fast verification of proof pillars P1, P2 with exact Fraction arithmetic.
P2: (d_a+d_b)E_n(u) = u E_n(u) (polynomial differentiation, exact).
P1: [z^E]g_{k,n} = P^E(d_a,d_b) applied before H-eval; checked for (n,k)=(2,4), E={2,2},{3},{4},{2,3}.
Uses step_02's machinery (E-series in (y,A,L), DB in directional cumulants) plus zeta-free DH.
"""
import sys, json
from fractions import Fraction as Fr
from math import comb, factorial
sys.path.insert(0, '.')
from step_02 import E_series, make_DB, poch_ratio, mono_name, H, const, padd, pmul, pscale, zeta_var, ZV

# ---- zeta-free DH(i,l,pj) at (n+j, n): Bell in kap^H (zeta parts removed) ----
def psiH(m, t):
    if m == 1:
        return const(H(1, t-1))
    f = Fr((-1)**m)
    for i in range(1, m): f *= i
    return const(-f*H(m, t-1))   # only the harmonic part

def make_DH(nv):
    def kapH(r, s, pj):
        p = nv+pj; q = nv; m = r+s
        out = pmul(const(Fr((-(nv-1))**r * (-1)**s)), psiH(m, p))
        out = padd(out, pmul(const(Fr(nv**r)), psiH(m, q)))
        if s == 0:
            out = padd(out, pscale(psiH(m, p+q), -1))
        return out
    cache = {}
    def DH(i, l, pj):
        key = (i, l, pj)
        if key in cache: return cache[key]
        if i == 0 and l == 0:
            res = const(1)
        elif i > 0:
            res = {}
            for r in range(i):
                for s in range(l+1):
                    T = pmul(const(Fr(comb(i-1, r)*comb(l, s))), pmul(kapH(r+1, s, pj), DH(i-1-r, l-s, pj)))
                    res = padd(res, T)
        else:
            res = {}
            for s in range(l):
                T = pmul(const(Fr(comb(l-1, s))), pmul(kapH(0, s+1, pj), DH(0, l-1-s, pj)))
                res = padd(res, T)
        cache[key] = res
        return res
    return DH

def Heval(nv, P, DH):
    # P: poly in (y,A,L) dict -> expand A^i L^l into a^r b^s and H-evaluate:
    # A = -(n-1)a + n b ; L = b - a
    tot = {}
    for (j, iA, iL), cv in P.items():
        # expand ((-(n-1))a + n b)^iA * (b-a)^iL into monomials a^r b^s
        for r in range(iA+1):
            cA = comb(iA, r) * (-(nv-1))**r * nv**(iA-r)   # a^r b^(iA-r)
            for s2 in range(iL+1):
                cL = comb(iL, s2) * (-1)**s2                 # a^s2 b^(iL-s2)  from (b-a)
                ra = r + s2
                sb = (iA-r) + (iL-s2)
                T = pmul(const(cv * poch_ratio(nv, j) * cA * cL), DH(ra, sb, j))
                tot = padd(tot, T)
    return tot

# ---- P2 ----
print('--- P2: (d_a+d_b)E_n = u E_n ---')
for nv in (1, 2, 3):
    Ek = E_series(nv, 6)
    ok = True
    for k in range(7):
        # (d_a+d_b) on Ek[k] should equal Ek[k-1] (coefficient of u^{k-1} in uE)
        D = {}
        for (j, iA, iL), v in Ek[k].items():
            # d_a + d_b = d_A  (since A = -(n-1)a+nb => d_aA+d_bA = 1; L = b-a => d_aL+d_bL = 0)
            if iA > 0:
                key = (j, iA-1, iL)
                D[key] = D.get(key, 0) + v*iA
        lhs = D
        rhs = Ek[k-1] if k >= 1 else {}
        if {kk: vv for kk, vv in padd(pscale(lhs, 1), pscale(rhs, -1)).items()}:
            ok = False
    print(f'n={nv}: (d_a+d_b)E = uE holds to order u^6: {ok}')

# ---- P1 ----
print()
print('--- P1: operator formula for (n,k)=(2,4) ---')
nv, k = 2, 4
Ek = E_series(nv, k)
DH = make_DH(nv)
# operators (in (a,b)); translate to (A,L): d_a = -(n-1) d_A - d_L ; d_b = n d_A + d_L
def apply_op(P, dA, dL):  # d_A^dA d_L^dL on poly
    R = {}
    for (j, iA, iL), v in P.items():
        if iA >= dA and iL >= dL:
            key = (j, iA-dA, iL-dL)
            fac = factorial(iA)//factorial(iA-dA) * factorial(iL)//factorial(iL-dL)
            R[key] = R.get(key, 0) + v*fac
    return R

def op_on_ab(P, ra, sb):  # d_a^ra d_b^sb : expand into (d_A, d_L)
    tot = {}
    for x in range(ra+1):
        for y_ in range(sb+1):
            # d_a^ra = sum_x C(ra,x) (-(n-1))^{ra-x} d_A^{ra-x} (-d_L)^x ; d_b^sb = sum C(sb,y) n^{sb-y} d_A^{sb-y} d_L^y
            pass
    # do it by composing monomial operators
    ops = {}
    for x in range(ra+1):
        for y2 in range(sb+1):
            cA = comb(ra, x) * (-(nv-1))**(ra-x) * (-1)**x
            cB = comb(sb, y2) * nv**(sb-y2)
            key = (ra-x + sb-y2, x + y2)
            ops[key] = ops.get(key, 0) + cA*cB
    tot = {}
    for (dA, dL), cv in ops.items():
        tot = padd(tot, pscale(apply_op(P, dA, dL), cv))
    return tot

def square(P, w):  # (d_a+d_b)^w - d_a^w - d_b^w = d_A^w - d_a^w - d_b^w  applied to P
    tot = apply_op(P, w, 0)  # d_A^w
    tot = padd(tot, pscale(op_on_ab(P, w, 0), -1))
    tot = padd(tot, pscale(op_on_ab(P, 0, w), -1))
    return tot

P4 = Ek[k]
t = json.load(open('mu_exact_N150_K6.json'))
pref = Fr(-nv*(2*nv-1))
def Cw(w): return Fr((-1)**(w+1)*factorial(w-1))
# E={2,2}: (1/2!)(C_2/2!)^2 square_2^2 ; E={3}: (C_3/3!) square_3 ; E={4}: (C_4/4!) square_4
tests = {
 'z2^2': (pscale(square(square(P4, 2), 2), (Cw(2)/2)**2 / 2), 'z2^2'),
 'z3':   (pscale(square(P4, 3), Cw(3)/Fr(6)), 'z3'),
 'z4':   (pscale(square(P4, 4), Cw(4)/Fr(24)), 'z4'),
 'z2z3': (pscale(square(square(P4, 2), 3), (Cw(2)/2)*(Cw(3)/Fr(6))), 'z2z3'),
}
for name, (Q, mon) in tests.items():
    val = Heval(nv, Q, DH)
    got = val.get(tuple(0 for _ in range(ZV)), Fr(0))  # zeta-free result (rational)
    exp = pref * Fr(t[str(k)][str(nv)][mon])
    print(f'[{mon}]g_{{4,2}}: operator formula = {got}  expected = {exp}  match: {got == exp}')
'''
open('step_05.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_05.py'],
                   capture_output=True, text=True, timeout=1200)
print(r.stdout)
print('ERR', r.stderr[-2000:])
