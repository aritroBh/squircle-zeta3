"""prove_E3b78.py -- PROOF of diagonal blindness (E3) at depths delta = 7, 8, all p.
Identical structure to TOW3 prove_E3b56.py / TOW2 prove_E3b34.py (proved there at delta<=6);
e-family is tow4lib.forced (FORCED family, PROVED), rho_pp(delta) = a_delta(2p,p) (PROVED Bell),
tau extended to n = 8 (tow4lib, verified exactly). e00_pp is anti-real: Dc = -conj(Delta).
The feed of Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s) (delta' <= delta)
into the child (p,p-1) at depth delta (DL2-direct + fold) vanishes identically.
Method: cached cancelled fractions + exact lcm-numerator (polynomial in p over Q(i) == 0).
Usage: python prove_E3b78.py DELTA
"""
import sympy as sp, time, sys
import tow4lib as T4

p = sp.Symbol('p')
i = sp.I

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    return T4.forced(delta, s, dd, mm, conj)

def tau1(n, conj=False):
    return (-1)**n * T4.tau2f(n) if conj else T4.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T4.tau1f(n) if conj else T4.tau2f(n)

rB0d = 2*i*(4*p-1)/(p+1)

_ehat = {}
def ehat(delta, s, conj=False):
    key = (delta, s, conj)
    if key in _ehat: return _ehat[key]
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * (rB0d if not conj else rB0d) / (4*p-dpp-1)
    out = sp.cancel(sp.expand((p+1)*tot))
    _ehat[key] = out
    return out

def rho_pp(delta, conj=False):
    return T4.forced(delta, 0, 2*p, p, conj)

_D = {}
def Delta(delta, conj=False):
    key = (delta, conj)
    if key in _D: return _D[key]
    D = rho_pp(delta, conj) - sum( (-1)**s * ehat(delta, s, conj) for s in range(0, delta//2+1) )
    out = sp.cancel(D)
    _D[key] = out
    return out

_Dc = {}
def Dc(delta):
    if delta in _Dc: return _Dc[delta]
    out = sp.cancel(-Delta(delta, conj=True))
    _Dc[delta] = out
    return out

_g = {}
def gfun(t):
    if t in _g: return _g[t]
    out = sp.cancel(sum( -Dc(dp)*tau1(t-dp, conj=True)/(4*p-3-dp) for dp in range(0, t+1) ))
    _g[t] = out
    return out

def feed_terms(delta):
    terms = []
    for dp in range(0, delta+1):
        terms.append(p*tau2(delta-dp)*Delta(dp)/(4*p-3-dp))
    for t in range(0, delta+1):
        c = (-i)**(-4-t) * i**(-6-delta)
        f = sp.Integer(1)
        for u in range(delta-t): f *= (4*p-6-t-u)
        binom = f/sp.factorial(delta-t)
        terms.append(-p * gfun(t) * c * binom)
    return terms

def monic(f):
    f = sp.expand(f)
    if f == 0: return f, sp.Integer(1)
    lc = sp.LC(sp.Poly(f, p))
    return sp.expand(f/lc), lc

def check_zero(exprs):
    rats = []
    for e in exprs:
        f = sp.cancel(e)
        rats.append(sp.fraction(f))
    fac_pow = {}; deninfo = []
    for num, den in rats:
        content, flist = sp.factor_list(den)
        cur = {}
        for f, e in flist:
            fm, lc = monic(f)
            content *= lc**e
            if fm == 1: continue
            fac_pow[fm] = max(fac_pow.get(fm, 0), e)
            cur[fm] = e
        deninfo.append((content, cur))
    DEN = sp.Integer(1)
    for f, e in fac_pow.items(): DEN *= f**e
    DEN = sp.expand(DEN)
    total = sp.Integer(0)
    for (num, den), (content, cur) in zip(rats, deninfo):
        quot = DEN
        for f, e in cur.items():
            for _ in range(e):
                q, r = sp.div(quot, f, p); assert r == 0; quot = q
        total += num * quot / content
    return sp.expand(total) == 0

if __name__ == '__main__':
    delta = int(sys.argv[1])
    t0 = time.time()
    terms = feed_terms(delta)
    print(f'  built {len(terms)} feed terms ({round(time.time()-t0,1)}s)', flush=True)
    is0 = check_zero(terms)
    print(f'E3 blindness delta={delta}: exact numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0
