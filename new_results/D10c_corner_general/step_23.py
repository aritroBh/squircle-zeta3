"""step_23.py -- debug: which product coefficients fail recognition in Q+Qc at 5pi/6?"""
import mpmath as mp
from ir_utils import find_relations
from fractions import Fraction
mp.mp.dps = 220
log2 = mp.log(2); z2 = mp.pi**2/6; ipi = mp.pi*1j
TOL = mp.mpf(10)**(-180)

def prod_coeffs(a):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    r1, r2, t1 = 2*c - 1, -1 - 2*c, 1 + 2*s
    one, mone = mp.mpf(1), mp.mpf(-1)
    poles = [one, mone, ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    bres = [mp.mpf(-1), mp.mpf(1), mp.mpf(-1), 2*c, -2*c]
    W = {}
    for rho, ar, br in zip(poles, ares, bres):
        for sig, sg in [(one, ar), (mone, ar), (ca, -(ar+br)), (r1, br), (r2, br)]:
            W[(rho, sig)] = W.get((rho, sig), 0) + sg
    logw = {rho: -(ar + br) for rho, ar, br in zip(poles, ares, bres)}
    prods = []
    def canon_prods(x, cf):
        if x > mp.mpf('0.5'):
            prods.append((-cf, 'logx', 'log1-x'))
        else:
            pass
    def reduce_li2(x, w):
        if abs(x) < TOL: return
        if abs(x - 1) < TOL: return
        if 0 < x < 1: canon_prods(x, -w)
        elif x < 0:
            canon_prods(x/(x-1), +w); prods.append((+w/2, 'log1-x2', 'log1-x2'))
        else:
            canon_prods(1/x, +w); prods.append((+w/2, 'logx', 'logx'))
    def Cfun(rho, sig):
        return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi
    for (rho, sig), w in W.items():
        if w == 0: continue
        if rho == sig:
            prods.append((w/2, 'lt1', 'lt1'))
            if rho != one: prods.append((-w/2, 'l1', 'l1'))
            continue
        reduce_li2((t1 - rho)/(sig - rho), w)
        prods.append((w, 'reC', 'lt1'))
        if rho != one:
            reduce_li2((1 - rho)/(sig - rho), -w)
            prods.append((-w, 'reC', 'l1'))
    for rho, wr in logw.items():
        if wr == 0: continue
        prods.append((wr, 'log2', 'lt1'))
        if rho != one: prods.append((-wr, 'log2', 'l1'))
    return prods, c

a = 5*mp.pi/6
prods, c = prod_coeffs(a)
def rec(x, c):
    rels = find_relations([x, mp.mpf(1), c], maxcoeff=10**6, resid_decades=25, max_report=1)
    if not rels: return None
    row, mc, res = rels[0]
    return (Fraction(-row[1], row[0]), Fraction(-row[2], row[0])), mc, abs(res)
seen = {}
for cf, la, lb in prods:
    k = mp.nstr(cf, 60)
    seen.setdefault(k, []).append((la, lb))
for k, v in sorted(seen.items(), key=lambda kv: mp.mpf(kv[0])):
    x = mp.mpf(k)
    r = rec(x, c)
    print(f"  {mp.nstr(x, 20):>24} -> {r[0] if r else 'FAIL'}   mc={r[1] if r else '-'}, resid={mp.nstr(r[2],2) if r else '-'}, e.g. {v[0]}")
