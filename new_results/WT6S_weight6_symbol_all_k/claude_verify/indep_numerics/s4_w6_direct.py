# Task (4): W6(k) directly from the full kernel pickles (degree-5 keys (5,0),(4,1),(3,2)):
#   W6(k) = int_b^a  sum_keys log^i t log^j(1-t) * sum_poles coeff_pole dv/(v-pole)
# compared with sum_m C(k-6,m-1) E_{5,m}.  Methods A (v-split, exp-substitution), B (t-variable),
# C (plain tanh-sinh in v, coarse).  The coefficients come straight from the pickles; the quadrature
# integrand is the summed kernel (NOT assembled from the 12 separate N's).
import sys, json, pickle
from fractions import Fraction
from math import comb
import mpmath as mp

DPS = 60
mp.mp.dps = DPS + 25
s2 = mp.sqrt(2)
a = 1 + s2
b = s2 - 1
KDIR = '<project>/research/08_theory/WTS_weight_k_structure/theory/'
KEYS = [(5, 0), (4, 1), (3, 2)]
POLE_TUP = {(0, 0): mp.mpf(0), (1, 1): a, (-1, -1): -a, (1, -1): 1 - s2}   # r0 + r1 sqrt2

def fr(x):
    x = Fraction(x)
    return mp.mpf(x.numerator) / x.denominator

def load_kernel(k):
    K = pickle.load(open(KDIR + f'kern_{k}.pkl', 'rb'))
    coef = {}   # key -> list of (pole value, coefficient mp number, is_a)
    for key in KEYS:
        d = K[key]
        assert not d['ipos'] and not d['ineg'], (k, key)
        assert all(Fraction(x) == 0 for x in d['const'])
        lst = []
        for (r0, r1), (c0, c1, c2, c3) in d['real'].items():
            assert Fraction(c2) == 0 and Fraction(c3) == 0, 'imaginary part present'
            r = fr(r0) + fr(r1) * s2
            c = fr(c0) + fr(c1) * s2
            lst.append((r, c, (Fraction(r0), Fraction(r1)) == (1, 1)))
        coef[key] = lst
    return coef

PTS = [mp.mpf(x) for x in (0, 1, 3, 8, 20, 50, 120, 250)]
SMAX = mp.mpf(30) * DPS / 10 + 80
PTS = [x for x in PTS if x < SMAX] + [SMAX]
LN2 = mp.log(2)
PTSB = [LN2] + PTS[1:]

# ---- Method A ----
def kernelA_left(coef, s):
    d = (1 - b) * mp.exp(-s)
    v = b + d
    c = (v + a) / (4 * v)
    logt = (mp.log(1 - b) - s) + mp.log(c)
    log1mt = mp.log1p(-(d * c))
    tot = mp.mpf(0)
    for (i, j) in KEYS:
        S = sum(cf / (v - r) for r, cf, isa in coef[(i, j)])
        tot += logt ** i * log1mt ** j * S
    return tot * d

def kernelA_right(coef, s):
    d = s2 * mp.exp(-s)
    v = a - d
    c = (v + b) / (4 * v)
    omt = d * c
    log1mt = (mp.log(s2) - s) + mp.log(c)
    logt = mp.log1p(-omt)
    tot = mp.mpf(0)
    for (i, j) in KEYS:
        S = mp.mpf(0)
        for r, cf, isa in coef[(i, j)]:
            S += (-cf) if isa else cf * d / (v - r)       # d/(v-a) = -1
        tot += logt ** i * log1mt ** j * S
    return tot

def W6_A(coef):
    L = mp.quad(lambda s: kernelA_left(coef, s), PTS, maxdegree=12)
    R = mp.quad(lambda s: kernelA_right(coef, s), PTS, maxdegree=12)
    return L + R

# ---- Method B ----
def vt(t):
    u = 2 * t - 1
    return u + mp.sqrt(u * u + 1)

def kernelB_left(coef, s):
    t = mp.exp(-s)
    v = vt(t)
    log1mt = mp.log1p(-t)
    jac = 4 * v * v / (v * v + 1) * t
    tot = mp.mpf(0)
    for (i, j) in KEYS:
        S = sum(cf / (v - r) for r, cf, isa in coef[(i, j)])
        tot += (-s) ** i * log1mt ** j * S
    return tot * jac

def kernelB_right(coef, s):
    u = mp.exp(-s)
    t = 1 - u
    v = vt(t)
    logt = mp.log1p(-u)
    jac = 4 * v * v / (v * v + 1) * u
    wa = -v * (v + b) / (v * v + 1)         # = w_a(t) * u
    tot = mp.mpf(0)
    for (i, j) in KEYS:
        S = mp.mpf(0)
        for r, cf, isa in coef[(i, j)]:
            S += cf * wa if isa else cf * jac / (v - r)
        tot += logt ** i * (-s) ** j * S
    return tot

def W6_B(coef):
    L = mp.quad(lambda s: kernelB_left(coef, s), PTSB, maxdegree=12)
    R = mp.quad(lambda s: kernelB_right(coef, s), PTSB, maxdegree=12)
    return L + R

# ---- Method C (coarse, plain v-quadrature) ----
def W6_C(coef, cdps=30):
    old = mp.mp.dps
    mp.mp.dps = cdps + 10
    a_, b_ = 1 + mp.sqrt(2), mp.sqrt(2) - 1
    def f(v):
        with mp.workdps(200):
            t = (v - b_) * (v + a_) / (4 * v)
            omt = (a_ - v) * (v + b_) / (4 * v)
            lt, l1 = mp.log(t), mp.log(omt)
            tot = mp.mpf(0)
            for (i, j) in KEYS:
                S = sum(cf / (v - r) for r, cf, isa in coef[(i, j)])
                tot += lt ** i * l1 ** j * S
            return tot
    val = mp.quad(f, [b_, mp.mpf(1), a_], maxdegree=10)
    mp.mp.dps = old
    return val

# ---- compare ----
Emine = json.load(open('E_values_mine.json'))['E_A']
Eref = json.load(open('<project>/research/08_theory/WT6S_weight6_symbol_all_k/e5m_values.json'))['values']
Em = {m: mp.mpf(Emine[str(m)]) for m in range(1, 6)}
Er = {m: mp.mpf(Eref[f'E5,{m}']) for m in range(1, 6)}

def Wformula(k, E):
    return sum(comb(k - 6, m - 1) * E[m] for m in range(1, 6))

ks = [8, 15] + [k for k in range(6, 19) if k not in (8, 15)]
for k in ks:
    coef = load_kernel(k)
    wa_ = W6_A(coef)
    wb_ = W6_B(coef)
    fm = Wformula(k, Em)
    fr_ = Wformula(k, Er)
    line = (f'k={k:2d}  W6(A)={mp.nstr(wa_, 40)}  |A-B|={mp.nstr(abs(wa_-wb_), 3)}  '
            f'|A - sum C(k-6,m-1)E_mine|={mp.nstr(abs(wa_-fm), 3)}  |A - sum..E_json|={mp.nstr(abs(wa_-fr_), 3)}')
    print(line, flush=True)
    if k in (8, 15):
        wc = W6_C(coef)
        mp.mp.dps = DPS + 25
        print(f'      k={k}: coarse C (plain tanh-sinh in v, 30 digit target): W6(C)={mp.nstr(wc, 28)}  |C-A|={mp.nstr(abs(wc-wa_), 3)}')
        print(f'      k={k}: sum_m C(k-6,m-1)E5,m (mine) = {mp.nstr(fm, 40)}')
        print(f'      k={k}: binomials C(k-6,m-1), m=1..5:', [comb(k - 6, m - 1) for m in range(1, 6)])
        print(f'      k={k}: relative size: |W6|={mp.nstr(abs(wa_), 6)}, magnitudes of terms:',
              [mp.nstr(comb(k - 6, m - 1) * Em[m], 8) for m in range(1, 6)])
