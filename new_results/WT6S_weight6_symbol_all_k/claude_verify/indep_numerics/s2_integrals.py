# Task (2): the 12 integrals N_{(i,j),r} = int_b^a log^i t log^j(1-t) dv/(v-r),  t=(v^2+2v-1)/(4v)
# Two independent quadratures (A: v-variable, split at v=1, exponential substitution at each log endpoint;
# B: t-variable, t = e^{-s} / 1-t = e^{-s}, Jacobian weights),  + coarse C: plain tanh-sinh in v.
import sys, json, time
from fractions import Fraction
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 25            # working precision (guard digits)
OUT = sys.argv[2] if len(sys.argv) > 2 else 'n_values.json'

s2 = mp.sqrt(2)
a = 1 + s2
b = s2 - 1
POLES = {'0': mp.mpf(0), 'a': a, '-a': -a, '-b': -b}      # -b = 1 - sqrt2
KEYS = [(5, 0), (4, 1), (3, 2)]
CELLS = [(kk, p) for kk in KEYS for p in POLES]            # 12 cells
SMAX = mp.mpf(30) * DPS / 10 + 80  # tail e^{-SMAX} s^5 negligible
PTS = [0, 1, 3, 8, 20, 50, 120, 250]
PTS = [mp.mpf(x) for x in PTS if x < SMAX] + [SMAX]
print('dps', DPS, 'working', mp.mp.dps, 'SMAX', mp.nstr(SMAX, 5), flush=True)

# ---------- sanity: the factorisations used ----------
def t_of_v(v): return (v * v + 2 * v - 1) / (4 * v)
for v in [b, mp.mpf('0.5'), 1, mp.mpf('1.7'), a]:
    t = t_of_v(v)
    assert abs(t - (v - b) * (v + a) / (4 * v)) < mp.mpf(10) ** (-DPS)
    assert abs(1 - t - (a - v) * (v + b) / (4 * v)) < mp.mpf(10) ** (-DPS)
print('t(b)=', mp.nstr(t_of_v(b), 5), ' t(1)=', mp.nstr(t_of_v(mp.mpf(1)), 5), ' t(a)=', mp.nstr(t_of_v(a), 12))

# ---------- Method A ----------
memoAL, memoAR = {}, {}

def A_left(s):
    """v in [b,1]: v-b = d = (1-b) e^{-s};  returns dict cell -> integrand wrt s (includes dv = d ds)"""
    if s in memoAL: return memoAL[s]
    d = (1 - b) * mp.exp(-s)
    v = b + d
    c = (v + a) / (4 * v)
    logt = (mp.log(1 - b) - s) + mp.log(c)
    t = d * c
    log1mt = mp.log1p(-t)
    res = {}
    for (i, j) in KEYS:
        base = logt ** i * log1mt ** j * d
        for p, r in POLES.items():
            res[((i, j), p)] = base / (v - r)
    memoAL[s] = res
    return res

def A_right(s):
    """v in [1,a]: a-v = d = sqrt2 e^{-s}; pole r=a: d/(v-a) = -1 exactly"""
    if s in memoAR: return memoAR[s]
    d = s2 * mp.exp(-s)
    v = a - d
    c = (v + b) / (4 * v)
    omt = d * c
    log1mt = (mp.log(s2) - s) + mp.log(c)
    logt = mp.log1p(-omt)
    res = {}
    for (i, j) in KEYS:
        base = logt ** i * log1mt ** j
        for p, r in POLES.items():
            if p == 'a':
                res[((i, j), p)] = -base
            else:
                res[((i, j), p)] = base * d / (v - r)
    memoAR[s] = res
    return res

def quadA(cell):
    fl = lambda s: A_left(s)[cell]
    fr = lambda s: A_right(s)[cell]
    L, el = mp.quad(fl, PTS, error=True, maxdegree=12)
    R, er = mp.quad(fr, PTS, error=True, maxdegree=12)
    return L + R, max(abs(el), abs(er))

# ---------- Method B ----------
memoBL, memoBR = {}, {}
LN2 = mp.log(2)
PTSB = [LN2] + [x for x in PTS if x > LN2]

def vt(t):
    u = 2 * t - 1
    return u + mp.sqrt(u * u + 1)

def B_left(s):
    """t = e^{-s}, s in [ln2, inf): int F(t) dt = int F(e^{-s}) e^{-s} ds"""
    if s in memoBL: return memoBL[s]
    t = mp.exp(-s)
    v = vt(t)
    log1mt = mp.log1p(-t)
    jac = 4 * v * v / (v * v + 1) * t
    res = {}
    for (i, j) in KEYS:
        base = (-s) ** i * log1mt ** j * jac
        for p, r in POLES.items():
            res[((i, j), p)] = base / (v - r)
    memoBL[s] = res
    return res

def B_right(s):
    """u = 1-t = e^{-s}; for pole a: w_a * u = -v (v+b)/(v^2+1)  (since 1/(v-a) = -(v+b)/(4 v (1-t)))"""
    if s in memoBR: return memoBR[s]
    u = mp.exp(-s)
    t = 1 - u
    v = vt(t)
    logt = mp.log1p(-u)
    jac = 4 * v * v / (v * v + 1) * u
    res = {}
    for (i, j) in KEYS:
        base = logt ** i * (-s) ** j
        for p, r in POLES.items():
            if p == 'a':
                res[((i, j), p)] = base * (-v * (v + b) / (v * v + 1))
            else:
                res[((i, j), p)] = base * jac / (v - r)
    memoBR[s] = res
    return res

def quadB(cell):
    fl = lambda s: B_left(s)[cell]
    fr = lambda s: B_right(s)[cell]
    L, el = mp.quad(fl, PTSB, error=True, maxdegree=12)
    R, er = mp.quad(fr, PTSB, error=True, maxdegree=12)
    return L + R, max(abs(el), abs(er))

# ---------- run ----------
NA, NB = {}, {}
t0 = time.time()
for cell in CELLS:
    NA[cell], eA = quadA(cell)
    NB[cell], eB = quadB(cell)
    print('N', cell, 'A=', mp.nstr(NA[cell], 30), ' |A-B|=', mp.nstr(abs(NA[cell] - NB[cell]), 3),
          ' errA', mp.nstr(eA, 3), ' errB', mp.nstr(eB, 3), ' t=%.0fs' % (time.time() - t0), flush=True)

worst = max(abs(NA[c] - NB[c]) for c in CELLS)
print('worst |A-B| over 12 integrals:', mp.nstr(worst, 5))

# ---------- Method C (coarse): plain tanh-sinh in v, factorised log arguments, inner extra precision ----------
mp.mp.dps = 40
CDPS = 30
def quadC(cell):
    (i, j), p = cell
    with mp.workdps(CDPS + 10):
        a_, b_, s2_ = 1 + mp.sqrt(2), mp.sqrt(2) - 1, mp.sqrt(2)
        r = {'0': mp.mpf(0), 'a': a_, '-a': -a_, '-b': -b_}[p]
        def f(v):
            with mp.workdps(200):
                t = (v - b_) * (v + a_) / (4 * v)
                omt = (a_ - v) * (v + b_) / (4 * v)
                return mp.log(t) ** i * mp.log(omt) ** j / (v - r)
        return mp.quad(f, [b_, mp.mpf(1), a_], maxdegree=10)

if len(sys.argv) > 3 and sys.argv[3] == 'C':
    tC = time.time()
    worstC = 0
    for cell in CELLS:
        mp.mp.dps = DPS + 25
        cv = quadC(cell)
        mp.mp.dps = DPS + 25
        dC = abs(cv - NA[cell])
        worstC = max(worstC, dC)
        print('C', cell, mp.nstr(cv, 25), ' |C-A|=', mp.nstr(dC, 3), flush=True)
    print('worst |C-A| (coarse, %d digit target):' % CDPS, mp.nstr(worstC, 5), ' t=%.0fs' % (time.time() - tC))

mp.mp.dps = DPS + 25
json.dump({'dps': DPS, 'A': {f'{c[0][0]}{c[0][1]}|{c[1]}': mp.nstr(NA[c], DPS + 24) for c in CELLS},
           'B': {f'{c[0][0]}{c[0][1]}|{c[1]}': mp.nstr(NB[c], DPS + 24) for c in CELLS}},
          open(OUT, 'w'), indent=1)
print('wrote', OUT)
