"""step_22.py -- PART 2b: reduction to canonical reps + LLL recognition of all
coefficients in [1,c] at alpha=5pi/6 (deg 4) and pi/12 (deg 8). dps 220.
Verifies the reduced closed form against a FRESH corner_series M1 evaluation.
"""
import mpmath as mp
from corner_series import corner_coeffs
from ir_utils import find_relations
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
            if abs(e[0] - rep) < TOL:
                e[1] += cf; return
        reps.append([rep, cf])
    zeta2 = [mp.mpf(0)]; prods = []; Jim = [mp.mpf(0)]
    def canon(x, cf):
        if x > mp.mpf('0.5'):
            add_li2(1 - x, -cf); zeta2[0] += cf
            prods.append((-cf, mp.log(x), mp.log(1 - x)))
        else:
            add_li2(x, cf)
    def reduce_li2(x, w):
        if abs(x) < TOL: return
        if abs(x - 1) < TOL:
            zeta2[0] -= w; return
        if 0 < x < 1: canon(x, -w)
        elif x < 0:
            canon(x/(x - 1), +w); prods.append((+w/2, mp.log(1 - x), mp.log(1 - x)))
        else:
            canon(1/x, +w); zeta2[0] -= 2*w
            prods.append((+w/2, mp.log(x), mp.log(x)))
            Jim[0] += w*mp.log(x)   # track imag (D10b convention: cancels)
    def Cfun(rho, sig):
        return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi
    for (rho, sig), w in W.items():
        if w == 0: continue
        if rho == sig:
            prods.append((w/2, mp.log(t1 - rho), mp.log(t1 - rho)))
            if rho != one: prods.append((-w/2, mp.log(1 - rho), mp.log(1 - rho)))
            continue
        C = Cfun(rho, sig)
        reduce_li2((t1 - rho)/(sig - rho), w)
        prods.append((w, mp.re(C), mp.log(t1 - rho)))
        if mp.im(C) != 0: Jim[0] -= w*mp.log(t1 - rho)
        if rho != one:
            reduce_li2((1 - rho)/(sig - rho), -w)
            prods.append((-w, mp.re(C), mp.log(1 - rho)))
            if mp.im(C) != 0: Jim[0] += w*mp.log(1 - rho)
    for rho, wr in logw.items():
        if wr == 0: continue
        prods.append((wr, log2, mp.log(t1 - rho)))
        if rho != one: prods.append((-wr, log2, mp.log(1 - rho)))
    total = mp.fsum([cf*mp.polylog(2, rep) for rep, cf in reps]) + zeta2[0]*z2 \
          + mp.fsum([cf*la*lb for cf, la, lb in prods])
    return reps, zeta2[0], prods, total, Jim[0], c

def rec_in_Qc(x, c, maxcoeff=10**6):
    rels = find_relations([x, mp.mpf(1), c], maxcoeff=maxcoeff, resid_decades=25, max_report=1)
    if not rels: return None
    row, mc, res = rels[0]
    from fractions import Fraction
    return (Fraction(-row[1], row[0]), Fraction(-row[2], row[0])), mc, res

for nm, a in [("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12)]:
    reps, zc, prods, total, Jim, c = reduced_closed_form(a)
    ref = corner_coeffs(a, 1)[1]
    res = total - ref
    dig = float(-mp.log10(abs(res))) if res != 0 else 999.0
    print(f"=== alpha = {nm} ===")
    print(f"  reduced closed form vs fresh M1: resid {mp.nstr(res,2)} ({dig:.1f} dig), J_imag {mp.nstr(Jim,2)}")
    ok = True
    reps.sort(key=lambda e: e[0])
    print(f"  Li2 coefficients recognition in Q+Qc (c = {mp.nstr(c,12)}):")
    for rep, cf in reps:
        r = rec_in_Qc(cf, c)
        tag = "OK" if r else "FAIL"
        ok &= bool(r)
        print(f"    Li2({mp.nstr(rep,18)}): coeff {mp.nstr(cf,14)} -> {r[0] if r else '??'} [{tag}]")
    r = rec_in_Qc(zc, c)
    ok &= bool(r)
    print(f"  zeta2 coeff {mp.nstr(zc,14)} -> {r[0] if r else '??'} [{'OK' if r else 'FAIL'}]")
    # product coefficients: unique cfs
    pcfs = sorted(set(mp.nstr(cf, 60) for cf, _, _ in prods))
    bad = 0
    for cfs in pcfs:
        if rec_in_Qc(mp.mpf(cfs), c) is None: bad += 1
    print(f"  product coefficients: {len(pcfs)} distinct, not-in-Q+Qc: {bad}")
    print(f"  ALL COEFFICIENTS IN Q+Qc: {ok and bad == 0}")
