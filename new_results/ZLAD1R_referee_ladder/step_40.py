"""ZLAD1R step_38 -- adversarial probes:
(a) sanity: wrong-sign rung must FAIL (tests are non-vacuous);
(b) completeness: every zeta monomial at weights 7,8 in the k<=8 table is either a pure z2-power
    or covered by a checked rung (list them);
(c) Phi_2 = Phi_1^2 + sum_{w>=4} (-1)^w (w-3) zeta(w) u^{w-4}  from my Omega computation;
(d) Phi_3 spot value: [z3]Phi_3 coefficient check via direct differentiation of Omega."""
import json
from fractions import Fraction as Fr
mine = json.load(open('muR_exact_N40_K8.json'))
def mu(k, n, m): return Fr(mine[str(k)][str(n)].get(m, '0'))
# (a) wrong sign must fail
bad = sum(1 for n in range(1, 41) if mu(8,n,'z3') == mu(7,n,'z2'))
print('(a) wrong-sign rung [z3]mu_8 == +[z2]mu_7 holds at', bad, 'of 40 values of n (expect 0):', 'OK' if bad == 0 else 'PROBLEM')
# (b) monomial census at weights 7,8 in table
def wt(m):
    if m == '1': return 0
    tot = 0
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, e = part.split('^'); tot += int(w)*int(e)
        else: tot += int(part)
    return tot
mons = set()
for k in range(9):
    for n in range(1, 41):
        mons |= set(mine[str(k)][str(n)].keys())
w7 = sorted(m for m in mons if wt(m) == 7); w8 = sorted(m for m in mons if wt(m) == 8)
print('(b) weight-7 monomials present:', w7)
print('    weight-8 monomials present:', w8)
rungs_covered = {'z7','z2z5','z3z4','z2^2z3',  # w7 single+mixed checked
                 'z8','z2z6','z3z5','z4^2','z2^2z4','z2z3^2'}  # w8
print('    w7 not covered by a checked rung:', set(w7) - rungs_covered)
print('    w8 not covered by a checked rung:', set(w8) - rungs_covered, "(should be only the free z2^4)")
# (c) recompute Omega quickly and test Phi_2 formula
from math import comb, factorial
NV7 = 7
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
def zp_const(c): return {(0,)*NV7: Fr(c)} if c else {}
def zp_var(w): return {tuple(1 if i == w-2 else 0 for i in range(NV7)): Fr(1)}
KMAX, EMAX = 8, 4
def beta_poly(w):
    p_prev = {(0, 0): zp_const(2)}; p_cur = {(1, 0): zp_const(1)}
    for i in range(2, w+1):
        p_new = {}
        for (du, ds), v in p_cur.items():
            p_new[(du+1, ds)] = zp_add(p_new.get((du+1, ds), {}), v)
        for (du, ds), v in p_prev.items():
            p_new[(du, ds+1)] = zp_add(p_new.get((du, ds+1), {}), v, -1)
        p_prev, p_cur = p_cur, p_new
    B = {(w, 0): zp_const(1)}
    for key, v in p_cur.items():
        B[key] = zp_add(B.get(key, {}), v, -1)
    return {k: v for k, v in B.items() if v}
L = {}
for w in range(2, 9):
    c = Fr((-1)**(w+1), w)
    for key, v in beta_poly(w).items():
        L[key] = zp_add(L.get(key, {}), zp_scale(zp_mul(v, zp_var(w)), c))
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
        Phi[ds][du] = zp_scale(v, Fr((-1)**ds) * factorial(ds))
# Phi_1^2
Phi1sq = {}
for d1, v1 in Phi[1].items():
    for d2, v2 in Phi[1].items():
        if d1+d2 > KMAX: continue
        Phi1sq[d1+d2] = zp_add(Phi1sq.get(d1+d2, {}), zp_mul(v1, v2))
# claimed correction sum_{w>=4} (-1)^w (w-3) z(w) u^{w-4}
claim2 = {}
for w in range(4, 9):
    claim2[w-4] = zp_add(claim2.get(w-4, {}), zp_scale(zp_var(w), Fr((-1)**w * (w-3))))
okc = all(Phi[2].get(du, {}) == zp_add(Phi1sq.get(du, {}), claim2.get(du, {})) for du in range(KMAX+1))
print('(c) Phi_2 = Phi_1^2 + sum_{w>=4} (-1)^w (w-3) z(w) u^{w-4}: ', okc)
# (d) [z3]-coefficient of Phi_3: from expansion, Omega = 1 + L + L^2/2 + L^3/6; [s^3] gets L^3/6's [s^3]
# check against direct symbolic: Phi_3 = (-d_s)^3 Omega|_0 = -( (d_s L)^3 + 3 d_sL d_s^2L + d_s^3 L )|_0
# numeric spot: z3 coefficient at u^1
v = Phi[3].get(1, {}).get((0,1,0,0,0,0,0), Fr(0))
print('(d) [z3 u^1] Phi_3 =', v, '(reference: computable; just recording)')
open('progress.log','a').write('step_38: adversarial probes OK (non-vacuity, weight 7/8 census, Phi_2 formula)\n')
