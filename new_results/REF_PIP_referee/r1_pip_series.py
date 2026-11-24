"""Referee script 1 (independent code): pi_p from its definition, and the 1/p expansion.

pi_p = (2/p) int_0^1 (t^{1-p} + (1-t)^{1-p})^{1/p} dt         (paper's definition)
Own derivation (checked by hand in REFEREE.md):
  pi_p = 4*2^{-h} + R,  h = 1/p,
  R = 4 h^2/(1-h) int_0^inf exp(-s h^2/(1-h)) (1+exp(-s h/(1-h)))^{-1-h} ((1+e^{-s})^h - 1) ds.
Term-by-term expansion in h (sympy, exact in symbols s, L=log(1+e^{-s}), lam=log 2), then
int s^i L^j ds = J(j,i) evaluated by quadrature at high precision.
"""
import sys, time, pickle
import mpmath as mp
import sympy as sp

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11   # highest N_k
DPS = 170
mp.mp.dps = DPS
T0 = time.time()
def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)

# ---------- 1. direct values of pi_p ----------
def pi_raw(p):
    p = mp.mpf(p)
    f = lambda t: (t**(1 - p) + (1 - t)**(1 - p))**(1 / p)
    return 2 / p * mp.quad(f, [0, mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4, 1])

def R_sform(p):
    p = mp.mpf(p); h = 1 / p
    c = h / (1 - h)
    f = lambda s: mp.exp(-s * h * c) * (1 + mp.exp(-s * c)) ** (-1 - h) * mp.expm1(h * mp.log1p(mp.exp(-s)))
    pts = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    return 4 * h * c * mp.quad(f, pts)

def pi_split(p):
    return 4 * mp.power(2, -1 / mp.mpf(p)) + R_sform(p)

mp.mp.dps = 60
log('pi_2 raw - pi =', mp.nstr(pi_raw(2) - mp.pi, 3))
for p in [mp.mpf(3), mp.mpf(5), mp.mpf(10), mp.mpf(1.5)]:
    a, b = pi_raw(p), pi_split(p)
    log('p=%s raw=%s split-raw=%s' % (mp.nstr(p, 5), mp.nstr(a, 45), mp.nstr(b - a, 3)))
# duality
for p in [3, 5, 10]:
    p = mp.mpf(p); q = p / (p - 1)
    log('duality p=%s: pi_p - pi_q (raw both) = %s' % (p, mp.nstr(pi_raw(p) - pi_raw(q), 3)))
mp.mp.dps = DPS

# ---------- 2. symbolic term-by-term expansion ----------
h, s, L, lam = sp.symbols('h s L lam')
K = KMAX + 2  # need h^{KMAX+1}
x = sp.symbols('x')
# log cosh(x) Taylor
lc = sp.series(sp.log(sp.cosh(x)), x, 0, K + 2).removeO()
def trunc(e, n=K):
    e = sp.expand(e)
    return sum(e.coeff(h, i) * h**i for i in range(n + 1))
a = trunc(sum(h**n for n in range(1, K + 2)))           # h/(1-h)
u = trunc(s * a)
g = lam - u / 2 + trunc(lc.subs(x, u / 2))               # log(1+e^{-u})
E = trunc(-s * h * a - (1 + h) * g)                      # log w  (constant term -lam)
E1 = trunc(E + lam)                                      # zero constant term
def series_exp0(z):  # exp of series with zero constant term
    out, term = sp.Integer(1), sp.Integer(1)
    for n in range(1, K + 1):
        term = trunc(term * z / n)
        out += term
    return trunc(out)
w = trunc(series_exp0(E1) / 2)                           # times e^{-lam} = 1/2
br = trunc(sum((h * L)**j / sp.factorial(j) for j in range(1, K + 1)))
integrand = trunc(trunc(h * a) * trunc(w * br))          # R/4 = int integrand ds
log('symbolic series done')
Jsym = {}
def to_J(c):
    c = sp.expand(c)
    if c == 0: return sp.Integer(0)
    poly = sp.Poly(c, s, L)
    out = 0
    for (i, j), co in poly.terms():
        assert j >= 1
        Jsym[(j, i)] = sp.Symbol('J_%d_%d' % (j, i))
        out += co * Jsym[(j, i)]
    return sp.expand(out)
Ncoef = {}
for k in range(0, KMAX + 1):
    Ncoef[k] = sp.expand((-lam)**(k + 1) / sp.factorial(k + 1) + to_J(integrand.coeff(h, k + 1)))
for k in range(0, 6):
    log('N_%d =' % k, Ncoef[k])
pickle.dump({'N': Ncoef, 'Jkeys': sorted(Jsym)}, open('r1_Nsym.pkl', 'wb'))

# ---------- 3. J numerics ----------
JV = {}
for (j, i) in sorted(Jsym):
    f = lambda t: t**i * mp.log1p(mp.exp(-t))**j
    JV[(j, i)] = mp.quad(f, [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 400])
log('J numerics done: %d integrals' % len(JV))
# exact J_{1,k} = k! eta(k+2), J_{2,0} = zeta(3)/4
for k in range(0, 6):
    if (1, k) in JV:
        log('J_1_%d - k! eta(k+2) = %s' % (k, mp.nstr(JV[(1, k)] - mp.factorial(k) * mp.altzeta(k + 2), 3)))
log('J_2_0 - zeta(3)/4 =', mp.nstr(JV[(2, 0)] - mp.zeta(3) / 4, 3))
subs = {lam: mp.log(2)}
def Nnum(k):
    e = Ncoef[k]
    tot = mp.mpf(0)
    for term in sp.Add.make_args(e):
        c, rest = term.as_coeff_Mul()
        val = mp.mpf(sp.Rational(c).p) / sp.Rational(c).q
        for fac in sp.Mul.make_args(rest):
            b, ex = fac.as_base_exp()
            if b == lam: val *= mp.log(2)**int(ex)
            elif b == 1: pass
            else:
                nm = str(b); _, j, i = nm.split('_'); val *= JV[(int(j), int(i))]**int(ex)
        tot += val
    return tot
NV = [Nnum(k) for k in range(KMAX + 1)]
for k in range(KMAX + 1):
    log('N_%d = %s' % (k, mp.nstr(NV[k], 45)))

# ---------- 4. printed closed forms ----------
l2, z3, z5, pi = mp.log(2), mp.zeta(3), mp.zeta(5), mp.pi
CF = {
 0: -l2, 1: l2**2 / 2, 2: pi**2 / 24 - l2**3 / 6,
 3: l2**4 / 24 + pi**2 / 24 * (1 - l2) + z3 / 4,
 4: -l2**5 / 120 + pi**2 / 24 * (1 - l2 + l2**2 / 2) + z3 / 4 * (1 - l2) + pi**4 / 2880,
 5: l2**6 / 720 + pi**2 / 24 * (1 - l2 + l2**2 / 2 - l2**3 / 6) + z3 / 4 * (1 - l2 + l2**2 / 2)
    + pi**4 / 2880 * (1 - l2) - 7 * pi**4 / 1920 - z3 * pi**2 / 48 + z5 / 8,
}
for k in range(6):
    log('closed form N_%d: diff = %s' % (k, mp.nstr(NV[k] - CF[k], 3)))
# E(p) expansion printed in the proof sketch
log('E coeff p^-3: 4*(N2 + log^3 2/6) - pi^2/6 =', mp.nstr(4 * (NV[2] + l2**3 / 6) - pi**2 / 6, 3))
log('E coeff p^-4: 4*(N3 - log^4 2/24) - (pi^2/6+z3-pi^2 log2/6) =',
    mp.nstr(4 * (NV[3] - l2**4 / 24) - (pi**2 / 6 + z3 - pi**2 * l2 / 6), 3))
# notes' 40-digit table
TAB = {6: '-0.9326682134840535164505582815863254270674', 7: '-1.472571756898086547743541305319946712566',
       8: '0.1221389382135391530907648898935863780245', 9: '4.700293562406214840901019828506766032978'}
for k, v in TAB.items():
    if k <= KMAX:
        log('notes table N_%d diff = %s' % (k, mp.nstr(NV[k] - mp.mpf(v), 3)))

# ---------- 5. direct pi_p vs truncated series; remainder order ----------
mp.mp.dps = 170
for p in [10**3, 10**4, 10**6, 10**8, 10**15]:
    P = mp.mpf(p); hh = 1 / P
    val = pi_split(P) / 4
    for KK in [3, 5, KMAX]:
        ser = 1 + sum(NV[k] * hh**(k + 1) for k in range(KK + 1))
        r = val - ser
        nxt = NV[KK + 1] * hh**(KK + 2) if KK + 1 <= KMAX else None
        log('p=1e%d  K=%d  pi_p/4 - series = %s   next term = %s' % (
            len(str(p)) - 1, KK, mp.nstr(r, 6), mp.nstr(nxt, 6) if nxt is not None else '-'))
# remainder order after N_3: (pi_p/4 - series_3)*p^5 -> N_4
for p in [10**2, 10**3, 10**4, 10**5]:
    P = mp.mpf(p); hh = 1 / P
    r = pi_split(P) / 4 - (1 + sum(NV[k] * hh**(k + 1) for k in range(4)))
    log('p=%d  (pi_p/4 - trunc_N3) * p^5 = %s   (N_4 = %s)' % (p, mp.nstr(r * P**5, 12), mp.nstr(NV[4], 12)))
# N_4, N_5 extracted from direct values alone (no series beyond N_3 / N_4): p=1e8
P = mp.mpf(10)**15; hh = 1 / P; val = pi_split(P) / 4
ext4 = (val - (1 + sum(CF[k] * hh**(k + 1) for k in range(4)) + sum(NV[k] * hh**(k + 1) for k in range(5, KMAX + 1)))) / hh**5
log('direct at p=1e15: (pi/4 - closed N0..N3)/h^5 - N4closed = %s  (N5.. subtracted numerically)' % (mp.nstr(ext4 - CF[4], 4),))
ext5 = (val - (1 + sum(CF[k] * hh**(k + 1) for k in range(6)) + sum(NV[k] * hh**(k + 1) for k in range(6, KMAX + 1))) + CF[5]*hh**6) / hh**6
log('direct at p=1e15: (pi/4 - closed N0..N4)/h^6 - N5closed = %s  (N6.. subtracted; expected ~N_{KMAX+1} h^{KMAX-4})' % (mp.nstr(ext5 - CF[5], 4),))
pickle.dump({k: str(v) for k, v in enumerate(NV)}, open('r1_Nnum.pkl', 'wb'))
log('done')
