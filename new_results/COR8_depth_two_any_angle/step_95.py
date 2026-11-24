
import sys, pickle
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 100
from cor8_series import *
from hpl_series import Li22
def li22_series(sgn):
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
S22m = li22_series(-1)   # Li22(1/2, 1-c)
S22p = li22_series(+1)   # Li22(1/2, 1+c)
DTS = sadd(sscale(sshift(S22m,-1), 1), sscale(S22m, -1), sscale(sshift(S22p,-1), -1), sscale(S22p, -1))
def DT(cc):
    return ((1-cc)/cc)*Li22(mp.mpf('0.5'), 1-cc) - ((1+cc)/cc)*Li22(mp.mpf('0.5'), 1+cc)
def evalS(S, cc):
    lc = mp.log(cc)
    return sum(v * cc**k * lc**l for (k,l), v in S.items())
for cs in ['1e-3','1e-4','1e-5']:
    cc = mp.mpf(cs)
    print(f"c={cs}: DT series vs numeric: {mp.nstr(abs(evalS(DTS,cc)-DT(cc)),3)}")
pickle.dump({'DTS': {k: str(v) for k,v in DTS.items()}}, open('series_dt.pkl','wb'))
