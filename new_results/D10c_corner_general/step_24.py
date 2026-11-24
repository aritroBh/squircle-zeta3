"""step_24.py -- PART 2b fixed: full-precision product-coefficient recognition at 5pi/6 and pi/12."""
import mpmath as mp
from corner_series import corner_coeffs
from ir_utils import find_relations
from fractions import Fraction
mp.mp.dps = 220
log2 = mp.log(2); z2 = mp.pi**2/6; ipi = mp.pi*1j
TOL = mp.mpf(10)**(-180)

def reduced_closed_form(a):
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
    reps = []
    def add_li2(rep, cf):
        for e in reps:
            if abs(e[0] - rep) < TOL: e[1] += cf; return
        reps.append([rep, cf])
    zeta2 = [mp.mpf(0)]; prods = []
    def canon(x, cf):
        if x > mp.mpf('0.5'):
            add_li2(1 - x, -cf); zeta2[0] += cf
            prods.append((-cf, mp.log(x), mp.log(1 - x)))
        else: add_li2(x, cf)
    def reduce_li2(x, w):
        if abs(x) < TOL: return
        if abs(x - 1) < TOL: zeta2[0] -= w; return
        if 0 < x < 1: canon(x, -w)
        elif x < 0: canon(x/(x - 1), +w); prods.append((+w/2, mp.log(1 - x), mp.log(1 - x)))
        else:
            canon(1/x, +w); zeta2[0] -= 2*w; prods.append((+w/2, mp.log(x), mp.log(x)))
    def Cfun(rho, sig):
        return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi
    for (rho, sig), w in W.items():
        if w == 0: continue
        if rho == sig:
            prods.append((w/2, mp.log(t1 - rho), mp.log(t1 - rho)))
            if rho != one: prods.append((-w/2, mp.log(1 - rho), mp.log(1 - rho)))
            continue
        reduce_li2((t1 - rho)/(sig - rho), w)
        prods.append((w, mp.re(Cfun(rho, sig)), mp.log(t1 - rho)))
        if rho != one:
            reduce_li2((1 - rho)/(sig - rho), -w)
            prods.append((-w, mp.re(Cfun(rho, sig)), mp.log(1 - rho)))
    for rho, wr in logw.items():
        if wr == 0: continue
        prods.append((wr, log2, mp.log(t1 - rho)))
        if rho != one: prods.append((-wr, log2, mp.log(1 - rho)))
    total = mp.fsum([cf*mp.polylog(2, rep) for rep, cf in reps]) + zeta2[0]*z2 \
          + mp.fsum([cf*la*lb for cf, la, lb in prods])
    return reps, zeta2[0], prods, total, c

def rec(x, c):
    rels = find_relations([x, mp.mpf(1), c], maxcoeff=10**6, resid_decades=30, max_report=1)
    if not rels: return None
    row, mc, res = rels[0]
    return (Fraction(-row[1], row[0]), Fraction(-row[2], row[0])), mc

for nm, a in [("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12)]:
    reps, zc, prods, total, c = reduced_closed_form(a)
    ref = corner_coeffs(a, 1)[1]
    dig = float(-mp.log10(abs(total - ref)))
    # dedupe product coefficients at full precision
    uniq = []
    for cf, la, lb in prods:
        if not any(abs(cf - u) < TOL for u in uniq): uniq.append(cf)
    fails = [(cf, rec(cf, c)) for cf in uniq]
    nfail = sum(1 for cf, r in fails if r is None)
    print(f"=== {nm}: total vs fresh M1 = {dig:.1f} digits; distinct product coeffs = {len(uniq)}, failures = {nfail}")
    for cf, r in sorted(fails, key=lambda z: z[0]):
        print(f"    {mp.nstr(cf, 16):>20} -> {r[0] if r else 'FAIL'}")
