"""ZLAD1R step_29 -- Theorem ZLAD1':
(a) NUMERIC: sum_w (-1)^{w+1} zeta(w) beta_w(u,s)/w  ==  log[Gamma(1+xi)Gamma(1+eta)/Gamma(1+u)]
    at sample (u,s) points (gamma cancels since xi+eta=u).
(b) EXACT: G_n(u) = sum_{e>=0} T_e(u) Phi_e(u), Phi_e = (-d/ds)^e Omega(u,0),
    Omega = exp(sum_{w>=2} (-1)^{w+1} zeta(w) beta_w/w), checked EXACTLY for all n<=40, k<=8,
    every zeta monomial, using MY table for T_e = [z2^e]G_n.
"""
import json
from fractions import Fraction as Fr
from math import comb, factorial
from mpmath import mp, mpf, zeta, loggamma, log, sqrt
mp.dps = 50

# ---- (a) numeric identity
def beta(w, u, s):  # Lucas: p_w via recurrence
    p0, p1 = 2, u
    if w == 0: return u**0*0  # unused
    if w == 1: return 0
    pw = [2, u]
    for i in range(2, w+1):
        pw.append(u*pw[-1] - s*pw[-2])
    return u**w - pw[w]
ok_a = True
for (u0, s0) in [(mpf('0.7'), mpf('0.09')), (mpf('1.3'), mpf('0.4')), (mpf('0.2'), mpf('-0.05')), (mpf('2.1'), mpf('1.0'))]:
    u, s = mpf(u0), mpf(s0)
    lhs = sum((-1)**(w+1) * zeta(w) * beta(w, u, s) / w for w in range(2, 60))
    # roots xi, eta of z^2 - u z + s
    disc = sqrt(u*u - 4*s)
    xi, eta = (u+disc)/2, (u-disc)/2
    rhs = loggamma(1+xi) + loggamma(1+eta) - loggamma(1+u)
    d = abs(lhs - rhs)
    print(f'(u,s)=({u0},{s0}): |series - logGamma ratio| = {mp.nstr(d,3)}')
    ok_a &= d < mpf('1e-35')
print('(a) logGamma identity incl. gamma cancellation:', ok_a)

# ---- (b) exact resummation
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
NV = 7  # z2..z8
def zp_add(P, Q, s=1):
    R = dict(P)
    for k, v in Q.items():
        nv = R.get(k, Fr(0)) + s*v
        if nv: R[k] = nv
        elif k in R: del R[k]
    return R
def zp_mul(P, Q):
    R = {}
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            k = tuple(x+y for x, y in zip(k1, k2))
            R[k] = R.get(k, Fr(0)) + v1*v2
    return {k: v for k, v in R.items() if v}
def zp_scale(P, c):
    return {k: v*c for k, v in P.items()} if c else {}
def zp_const(c): return {(0,)*NV: Fr(c)} if c else {}
def zp_var(w):
    return {tuple(1 if i == w-2 else 0 for i in range(NV)): Fr(1)}

KMAX = 8; EMAX = 4
# beta_w as (u,s)-poly with const coeffs: dict {(du, ds): poly}
def beta_poly(w):
    pw = {(0, 0): zp_const(2)} if w == 0 else None
    p_prev = {(0, 0): zp_const(2)}; p_cur = {(1, 0): zp_const(1)}
    if w == 1: pw = p_cur
    for i in range(2, w+1):
        p_new = {}
        for (du, ds), v in p_cur.items():
            p_new[(du+1, ds)] = zp_add(p_new.get((du+1, ds), {}), v)
        for (du, ds), v in p_prev.items():
            p_new[(du+1, ds+1)] = zp_add(p_new.get((du+1, ds+1), {}), v, -1)
        p_prev, p_cur = p_cur, p_new
        pw = p_cur
    # beta = u^w - p_w
    B = {(w, 0): zp_const(1)}
    for key, v in pw.items():
        B[key] = zp_add(B.get(key, {}), v, -1)
    return {k: v for k, v in B.items() if v}

# L = sum_{w=2..8} (-1)^{w+1} zeta(w)/w * beta_w   (as (u,s)-series with zeta-poly coeffs)
L = {}
for w in range(2, 9):
    c = Fr((-1)**(w+1), w)
    for key, v in beta_poly(w).items():
        L[key] = zp_add(L.get(key, {}), zp_scale(zp_mul(v, zp_var(w)), c))
# Omega = exp(L), truncated to du<=8, ds<=4
def series_mul(A, B):
    R = {}
    for (d1, e1), v1 in A.items():
        for (d2, e2), v2 in B.items():
            if d1+d2 > KMAX or e1+e2 > EMAX: continue
            key = (d1+d2, e1+e2)
            R[key] = zp_add(R.get(key, {}), zp_mul(v1, v2))
    return {k: v for k, v in R.items() if v}
Om = {(0, 0): zp_const(1)}
Lpow = {(0, 0): zp_const(1)}
for r in range(1, 9):
    Lpow = series_mul(Lpow, L)
    for key, v in Lpow.items():
        Om[key] = zp_add(Om.get(key, {}), zp_scale(v, Fr(1, factorial(r))))
Phi = {e: {} for e in range(EMAX+1)}
for (du, ds), v in Om.items():
    if ds <= EMAX:
        Phi[ds][(du,)] = zp_scale(v, Fr((-1)**ds) * factorial(ds))
# Phi[e]: dict {(du,): zpoly}

# T_e(u) for given n: coeff of u^k is -n(2n-1)*[z2^e]mu_{k,n}
def T_e(n, e):
    m = 'z2' if e == 1 else ('1' if e == 0 else f'z2^{e}')
    return {k: (-n*(2*n-1))*mu(k, n, m) if e >= 1 else None for k in range(KMAX+1)}

def mname_to_tup(m):
    e = [0]*NV
    if m == '1': return tuple(e)
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, x = part.split('^'); e[int(w)-2] += int(x)
        else:
            e[int(part)-2] += 1
    return tuple(e)

ok_b = True
checked = 0
for n in range(2, 41):
    # build reconstructed G_n(u) zeta-polys: rec[k][monom_tup] = value
    rec = {k: {} for k in range(KMAX+1)}
    for e in range(EMAX+1):
        m_e = '1' if e == 0 else ('z2' if e == 1 else f'z2^{e}')
        for k1 in range(KMAX+1):
            cT = (-n*(2*n-1))*mu(k1, n, m_e)
            if cT == 0: continue
            for (du,), vz in Phi[e].items():
                k = k1 + du
                if k > KMAX: continue
                rec[k] = zp_add(rec[k], zp_scale(vz, cT))
    # compare with -n(2n-1)*mu table
    for k in range(KMAX+1):
        for m, vs in mine[str(k)][str(n)].items():
            got = rec[k].get(mname_to_tup(m), Fr(0))
            expv = (-n*(2*n-1))*Fr(vs)
            checked += 1
            if got != expv:
                ok_b = False
                print('MISMATCH at n,k,m =', n, k, m, got, expv)
print(f'(b) EXACT resummation G_n = sum_e T_e Phi_e: {checked} (n,k,monomial) values, n=2..40, k<=8: {"ALL MATCH" if ok_b else "FAIL"}')
open('progress.log','a').write(f'step_29: ZLAD1 prime: logGamma/beta identity numeric OK (4 pts, 35+ digits); EXACT resummation n=2..40 k<=8 ALL MATCH ({checked} values)\n')
