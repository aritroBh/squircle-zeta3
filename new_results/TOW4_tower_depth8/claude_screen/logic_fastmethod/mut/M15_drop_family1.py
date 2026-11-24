"""prove_adj78_fast.py -- adjacent Newton identities at depths 7,8, all p (symbolic).
FAST exact method with CACHED cancelled Ediag/ehat fractions + exact lcm-numerator check.
Identity lhs == sum(terms) proved as exact rational identity in p over Q(i): numerator == 0.
Usage: python prove_adj78_fast.py DELTA S
"""
import sys, time
import sympy as sp
import tow4lib as T4

p = sp.Symbol('p')
d = 2*p-1
i = sp.I

def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    return T4.forced(delta, s, dd, mm, conj)

def tau1(n, conj=False):
    return (-1)**n * T4.tau2f(n) if conj else T4.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T4.tau1f(n) if conj else T4.tau2f(n)

rA0 = i*(4*p-3)/(p+1)
rB0 = (4*p-1)*(4*p-3)/((p+1)*p)

_ehat = {}
def ehat(delta, s, conj=False):
    key = (delta, s, conj)
    if key in _ehat: return _ehat[key]
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * rB0 / (4*p-dpp-1)
    out = sp.cancel(sp.expand((p+1)*tot))
    _ehat[key] = out
    return out

_ed = {}
def Ediag(delta, s, conj=False):
    key = (delta, s, conj)
    if key in _ed: return _ed[key]
    tot = 0; t = 0
    while s+t <= delta//2:
        tot += (-1)**t * ehat(delta, s+t, conj)
        t += 1
    out = sp.cancel(sp.expand(tot))
    _ed[key] = out
    return out

def feed_terms(delta, s):
    terms = []
    for dp in range(0, delta+1):
        pass
    for dp in range(0, delta+1):
        terms.append(p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))
    for t in range(0, delta+1):
        mu = 4*p-4-t
        g = 0
        for dp in range(0, t+1):
            g -= Ediag(dp, s, conj=True)*tau1(t-dp, conj=True)/(4*p-3-dp)
        binom = ff(mu-2, delta-t)/sp.factorial(delta-t)
        terms.append(-p * sp.expand(g) * (-i)**(-4-t) * binom * i**(-6-delta))
    return terms

def monic(f):
    f = sp.expand(f)
    if f == 0: return f, sp.Integer(1)
    lc = sp.LC(sp.Poly(f, p))
    return sp.expand(f/lc), lc

def main(delta, s):
    t0 = time.time()
    lhs = sp.cancel(eR(delta, s, d, p) + eR(delta, s+1, d, p))
    terms = feed_terms(delta, s)
    print(f'  built lhs + {len(terms)} terms ({round(time.time()-t0,1)}s)', flush=True)
    exprs = [lhs] + [-t for t in terms]
    rats = []
    for e in exprs:
        f = sp.cancel(e)
        rats.append(sp.fraction(f))
    print(f'  combined to single fractions ({round(time.time()-t0,1)}s)', flush=True)
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
    total = sp.expand(total)
    is0 = (total == 0)
    print(f'adjacent ({delta},{s}) FAST+cached: exact numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
