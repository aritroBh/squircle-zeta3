
import sys, json, re, pickle
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 100
from fractions import Fraction as Fr
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle/'
R3 = json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
from hpl_series import hpl_coeffs
om_terms = []
for name, coef in R3:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    if len(w) == 4:
        om_terms.append((e, w, Fr(coef)))
# Omega series (rational)
OMS = {}
for e, w, q in om_terms:
    cf = hpl_coeffs(w, K+5)
    S = {(k, 0): mp.mpf(0) for k in range(K+1)}
    for k in range(1, K+1):
        S[(k,0)] = cf[k]
    S = sscale(S, mp.mpf(q.numerator)/q.denominator)
    OMS = sadd(OMS, sshift(S, e))
print("Omega series first coeffs:", [(k, mp.nstr(OMS.get((2*k,0),0), 8)) for k in range(6)])
# DT series: Li22(1/2, 1-c), Li22(1/2, 1+c)
def li22_series(sgn):
    # [c^k] = sgn^k sum_{n>m>=k} C(m,k) 2^{-n}/(n^2 m^2)
    out = {}
    from math import comb
    for k in range(0, K+1):
        tot = mp.mpf(0)
        n = max(k+1, 2)
        while True:
            inner = mp.mpf(0)
            for m in range(max(k,1), n):
                inner += mp.mpf(comb(m,k))/m**2
            term = mp.power(2, -n)/n**2 * inner
            tot += term
            if n > 60 and abs(term) < mp.mpf(10)**(-mp.mp.dps-10):
                break
            n += 1
        out[(k, 0)] = tot * (sgn**k if k else 1)
    return trim(out)
S22m = li22_series(1)
S22p = li22_series(-1)
# DT = ((1-c)/c) S22m - ((1+c)/c) S22p
DTS = sadd(sscale(sshift(S22m,-1), 1), sscale(S22m, -1), sscale(sshift(S22p,-1), -1), sscale(S22p, -1))
print("DT series first coeffs:", [(k, mp.nstr(DTS.get((k,0),0), 8)) for k in range(-1, 6)])
pickle.dump({'OMS': {k: str(v) for k,v in OMS.items()}, 'DTS': {k: str(v) for k,v in DTS.items()}}, open('series_om_dt.pkl','wb'))
# sanity vs numeric at c=1e-3
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import H_series, Li22
def Omega(cc):
    return sum(mp.mpf(q.numerator)/q.denominator * cc**e * H_series(w, cc) for e, w, q in om_terms)
def DT(cc):
    return ((1-cc)/cc)*Li22(mp.mpf('0.5'), 1-cc) - ((1+cc)/cc)*Li22(mp.mpf('0.5'), 1+cc)
def evalS(S, cc):
    lc = mp.log(cc)
    return sum(v * cc**k * lc**l for (k,l), v in S.items())
for cs in ['1e-3','1e-4']:
    cc = mp.mpf(cs)
    print(f"c={cs}: Omega series vs numeric: {mp.nstr(abs(evalS(OMS,cc)-Omega(cc)),3)}, DT: {mp.nstr(abs(evalS(DTS,cc)-DT(cc)),3)}")
