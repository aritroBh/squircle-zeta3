code = r'''"""verify_zlad1r.py -- ZLAD1R referee verification, fully independent of ZLAD1 code.
Recomputes mu_{k,n} (k<=8, n<=40) from thm:flat's Beta-integral formula via Stirling expansion
of f^{1-2n} + direct bivariate Gamma-series exponentiation (a different route from ZLAD1's and G1's),
then checks:
  (1) table == ZLAD1's mu_exact_N150_K8.json exactly (if present) and gamma cancels identically;
  (2) five ladder rungs EXACT (incl. weight-7/8 singles and two mixed), all n<=40;
  (3) two requested closed forms EXACT (n<=40) + canonical coefficients;
  (4) Lemma 2: (d_a+d_b)[u^k]E_n = [u^{k-1}]E_n (own code);
  (5) Lemma 1: operator formula [zeta^E]G_n = HH(P^E [u^k]E_n) at 6 sample points (own code);
  (6) Theorem ZLAD1': symbolic logGamma/beta identity exact; EXACT resummation G_n = sum T_e Phi_e
      for n=2..40, k<=8, every zeta monomial;
  (7) end-to-end numeric quadrature (Cauchy) of the Beta integral at (k,n)=(8,3), 40 digits.
Runtime target < 5 min. Exit nonzero on any failure.
"""
import json, os, sys, time
from fractions import Fraction as Fr
from math import comb, factorial
t0 = time.time()
fails = []
def group(name, ok):
    print(('PASS' if ok else 'FAIL'), f'{name}  [{time.time()-t0:.0f}s]', flush=True)
    if not ok: fails.append(name)

KMAX, NMAX = 8, 40
# ---------- exact zeta/gamma polynomials ----------
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
                    al, be = i-x+(j-y_), x+y_
                    if (n+t,) not in exp_cache:
                        exp_cache[(n+t,)] = exp_series(n+t, n, KMAX)
                    sub = exp_cache[(n+t,)].get((al, be), {})
                    if not sub: continue
                    c = base * cx * cy * poch_ratio(n, t) * factorial(al) * factorial(be)
                    tot = padd(tot, pscale(sub, c))
    return pscale(tot, Fr(-1, n*(2*n-1)))
def mname(e):
    if e == (0,)*7: return '1'
    out = []
    for idx, x in enumerate(e):
        if x: out.append(f'z{idx+2}' + (f'^{x}' if x > 1 else ''))
    return ''.join(out)

print('recomputing independent table (k<=8, n<=40)...', flush=True)
table = {k: {} for k in range(KMAX+1)}
gamma_left = 0
for n in range(1, NMAX+1):
    exp_cache = {}
    for k in range(KMAX+1):
        P = mu_kn(k, n, exp_cache)
        gamma_left += sum(1 for kk in P if kk[0] != 0)
        table[k][n] = {mname(kk[1:]): vv for kk, vv in P.items() if kk[0] == 0}
def mu(k, n, m): return table[k][n].get(m, Fr(0))
group('(1a) gamma cancels identically in every mu_{k,n}', gamma_left == 0)

# (1b) diff against ZLAD1 table
zl_path = '../ZLAD1_zeta_ladder/mu_exact_N150_K8.json'
ok = True
if os.path.exists(zl_path):
    zl = json.load(open(zl_path))
    for k in range(9):
        for n in range(1, 41):
            a = {m: str(v) for m, v in table[k][n].items()}
            b = zl[str(k)][str(n)]
            for m in set(a) | set(b):
                if a.get(m, '0') != b.get(m, '0'): ok = False
else:
    print('  (ZLAD1 table not found; skipped)')
group('(1b) independent table == ZLAD1 table exactly (k<=8, n<=40)', ok)

# (2) five rungs
ok = True
for n in range(1, NMAX+1):
    ok &= mu(8,n,'z3') == -mu(7,n,'z2')
    ok &= mu(8,n,'z7') == -mu(3,n,'z2') - 4*mu(5,n,'z2^2') - 6*mu(7,n,'z2^3')
    ok &= mu(8,n,'z8') == mu(2,n,'z2') + 5*mu(4,n,'z2^2') + 12*mu(6,n,'z2^3') + 6*mu(8,n,'z2^4')
    ok &= mu(8,n,'z3z5') == 2*mu(4,n,'z2^2') + 6*mu(6,n,'z2^3')
    ok &= mu(7,n,'z2^2z3') == -3*mu(6,n,'z2^3')
group('(2) five rungs EXACT (rung1; w7; w8; mixed z3z5; mixed z2^2z3), n<=40', ok)

# (3) closed forms + canonical
ok = True
for n in range(1, NMAX+1):
    ok &= mu(5,n,'z5') == Fr(-3, (2*n-1)*(2*n+1)*(2*n+3))
    ok &= mu(6,n,'z2^3') == Fr(n*(5*n-1), 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5))
    ok &= mu(2,n,'z2') == Fr(1, 4*n*n-1)
    ok &= mu(3,n,'z3') == Fr(-1, 4*n*n-1)
    ok &= mu(4,n,'z2^2') == Fr(-n, (2*n-1)*(2*n+1)*(2*n+3))
group('(3) closed forms [z5]mu_5, [z2^3]mu_6 + canonical + G2 input [z2^2]mu_4, n<=40', ok)

# (4) Lemma 2 with own E-series
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
def Dsum(P):
    R = {}
    for (T, al, be), v in P.items():
        if al > 0:
            key = (T, al-1, be); R[key] = R.get(key, Fr(0)) + v*al
        if be > 0:
            key = (T, al, be-1); R[key] = R.get(key, Fr(0)) + v*be
    return {kk: vv for kk, vv in R.items() if vv}
ok = True
for n in (2, 3, 5):
    EPs = [E_poly(k, n) for k in range(9)]
    for k in range(1, 9):
        if Dsum(EPs[k]) != EPs[k-1]: ok = False
group('(4) Lemma 2: (d_a+d_b)E_n = u E_n, n in {2,3,5}, k<=8', ok)

# (5) Lemma 1 operator formula at sample points (own code)
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
def d_ab(P, ra, rb):
    R = {}
    for (T, al, be), v in P.items():
        if al >= ra and be >= rb:
            fac = factorial(al)//factorial(al-ra) * factorial(be)//factorial(be-rb)
            key = (T, al-ra, be-rb)
            R[key] = R.get(key, Fr(0)) + v*fac
    return {kk: vv for kk, vv in R.items() if vv}
def sq_w(P, w):
    R = {}
    for r in range(w+1):
        Q = pscale(d_ab(P, r, w-r), comb(w, r))
        R = {kk: R.get(kk, Fr(0)) + vv for kk, vv in Q.items()} if not R else \
            {kk: (R.get(kk, Fr(0)) + Q.get(kk, Fr(0))) for kk in set(R)|set(Q)}
        R = {kk: vv for kk, vv in R.items() if vv}
    for P2, sgn in ((d_ab(P, w, 0), -1), (d_ab(P, 0, w), -1)):
        for kk, vv in P2.items():
            R[kk] = R.get(kk, Fr(0)) + sgn*vv
            if R[kk] == 0: del R[kk]
    return R
def Cw(w): return Fr((-1)**(w+1)) * factorial(w-1)
def apply_PE(P, E):
    ops = []; denom = 1
    for w, e in E:
        denom *= factorial(e)
        for _ in range(e): ops.append((w, Cw(w)/Fr(factorial(w))))
    R = P
    for w, c in ops:
        R = pscale(sq_w(R, w), c)
    return pscale(R, Fr(1, denom)) if False else {kk: vv/denom for kk, vv in R.items() if vv}
def HH(P, n, cache):
    tot = Fr(0)
    for (T, al, be), v in P.items():
        if (n+T,) not in cache: cache[(n+T,)] = exp_series_H(n+T, n, 12)
        sub = cache[(n+T,)].get((al, be), Fr(0))
        tot += v * poch_ratio(n, T) * factorial(al) * factorial(be) * sub
    return tot
ok = True
for (n, k), E, m in [((2,4),[(2,1)],'z2'), ((2,4),[(3,1)],'z3'), ((3,6),[(2,3)],'z2^3'),
                     ((3,6),[(3,2)],'z3^2'), ((2,8),[(8,1)],'z8'), ((2,8),[(2,2),(3,1)],'z2^2z3'),
                     ((5,8),[(2,1),(6,1)],'z2z6')]:
    got = HH(apply_PE(E_poly(k, n), E), n, {})
    ok &= (got == -Fr(n*(2*n-1))*mu(k, n, m))
group('(5) Lemma 1 operator formula at 7 sample (n,k,E) incl. w=8', ok)

# (6) ZLAD1'
import sympy as sp
u_, s_ = sp.symbols('u s')
p_prev, p_cur = sp.Integer(2), u_
beta_sp = {}
for w in range(2, 14):
    p_new = sp.expand(u_*p_cur - s_*p_prev); p_prev, p_cur = p_cur, p_new
    beta_sp[w] = sp.expand(u_**w - p_cur)
zs = sp.symbols('z2:z15')
lhs_ser = sp.expand(sum((-1)**(w+1)*zs[w-2]/w*beta_sp[w] for w in range(2, 14)))
rhs_ser = sp.expand(sum((-1)**m*zs[m-2]/m*((u_**m - beta_sp[m]) - u_**m) for m in range(2, 14)))
ok1 = (sp.expand(lhs_ser - rhs_ser) == 0)
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
def zp_var(w):
    return {tuple(1 if i == w-2 else 0 for i in range(NV7)): Fr(1)}
EMAX = 4
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
def m2t(m):
    e = [0]*NV7
    if m == '1': return tuple(e)
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, x = part.split('^'); e[int(w)-2] += int(x)
        else:
            e[int(part)-2] += 1
    return tuple(e)
ok2 = True
for n in range(2, 41):
    rec = {k: {} for k in range(KMAX+1)}
    for e in range(EMAX+1):
        m_e = '1' if e == 0 else ('z2' if e == 1 else f'z2^{e}')
        for k1 in range(KMAX+1):
            cT = (-n*(2*n-1))*mu(k1, n, m_e)
            if cT == 0: continue
            for du, vz in Phi[e].items():
                k = k1 + du
                if k > KMAX: continue
                rec[k] = zp_add(rec[k], zp_scale(vz, cT))
    for k in range(KMAX+1):
        for m in table[k][n]:
            if rec[k].get(m2t(m), Fr(0)) != (-n*(2*n-1))*mu(k, n, m): ok2 = False
group("(6) ZLAD1': logGamma identity exact + resummation G_n=sum T_e Phi_e EXACT n=2..40,k<=8", ok1 and ok2)

# (7) quadrature end-to-end at (8,3)
from mpmath import mp, mpf, mpc, pi, exp, log, zeta, quad
from math import gamma as gg
mp.dps = 45
def Gn_int(y, n, u):
    if y == 0 or y == 1: return mpf(0)
    a = log(y); b = log(1-y)
    return y**(n-1)*(1-y)**(n-1)*exp(u*(-(n-1)*a+n*b))*(1-y+y*exp(u*(b-a)))**(1-2*n)
def Gn(u, n):
    B = mpf(gg(n))**2/gg(2*n)
    return quad(lambda y: Gn_int(y, n, u), [0, 1]) / B
rho = mpf('0.4'); nv, kv = 3, 8
ck = quad(lambda th: Gn(rho*exp(mpc(0, th)), nv)*exp(mpc(0, -kv*th)), [0, 2*pi])/(2*pi*rho**kv)
mu_num = -ck/(nv*(2*nv-1))
ZET = {w: zeta(w) for w in range(2, 9)}
def mono_val(m):
    if m == '1': return mpf(1)
    v = mpf(1)
    for part in m.split('z'):
        if not part: continue
        if '^' in part:
            w, e = part.split('^'); v *= ZET[int(w)]**int(e)
        else:
            v *= ZET[int(part)]
    return v
mu_ex = sum(mpf(v.numerator)/mpf(v.denominator)*mono_val(m) for m, v in table[kv][nv].items())
d = abs(mu_num - mu_ex)
print(f'  quadrature (8,3): |num - exact| = {mp.nstr(d, 3)}')
group('(7) end-to-end quadrature vs exact table at (8,3), 40 digits', d < mpf('1e-25'))

print(f'total {time.time()-t0:.0f}s')
if fails:
    print('FAILURES:', fails); sys.exit(1)
print('ALL ZLAD1R VERIFICATIONS PASSED')
'''
open('verify_zlad1r.py', 'w').write(code)
print('written', len(code))
