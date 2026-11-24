
step10 = r'''"""step_10.py -- G2b: the numerical (3,2) cobracket at the hexagon, including the zeta(3)-part.

delta_word (M4c's formula, exact-checked at the symbol level in step_05) with exact Lambda^2
coordinates and numerical 3-letter integrals from J3val (closed-form J2 + own tanh-sinh;
tangential regularization for coincident letters, PV at interior poles; validated: from-0 vs
series 25/26 + all-zero word, shuffle identity 8.5e-66, J2 closed form vs series/quad 1e-102).
D_num(P) fitted over {z3, Li3(1/3), Li3(-1/2), weight-3 products} by LLL at two precisions:
Li3 coefficients must match the EXACT symbol values (chain validation); z3 coefficient is new.
"""
import sys, os, json, time
from fractions import Fraction as Fr
import mpmath as mp
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, G2B)
from g2b_dbg import s05
from math import gcd
from functools import reduce
from flint import fmpz_mat, fmpz

DPS = 110
mp.mp.dps = DPS + 20
t0 = time.time()
C0 = Fr(1, 2)

def log(msg):
    print(msg, flush=True)
    with open(G2B + "progress.log", "a") as f:
        f.write(msg + "\n")

# ---------- validated J3 evaluator ----------
def reli2(z):
    z = mp.mpf(z)
    if z == 0: return mp.mpf(0)
    if z == 1: return mp.pi**2 / 6
    if abs(z) <= mp.mpf('0.55'):
        tot = mp.mpf(0); zk = z; k = 1
        while True:
            t = zk / (k * k); tot += t
            if abs(t) < mp.mpf(10) ** (-mp.mp.dps + 10): break
            zk *= z; k += 1
        return tot
    if z > 1: return -reli2(1 / z) + mp.pi**2 / 3 - mp.log(z)**2 / 2
    if z > 0: return mp.pi**2 / 6 - mp.log(z) * mp.log(1 - z) - reli2(1 - z)
    w = z / (z - 1)
    return -mp.log(1 - z)**2 / 2 - reli2(w)

def Gfun(u, b, c):
    if u == c: return mp.mpf(0)
    if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + reli2((u - c) / (b - c))

def J2val(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

def ts(f, a, b, level=5):
    if a == b: return mp.mpf(0)
    sgn = 1
    if b < a: a, b = b, a; sgn = -1
    h = mp.mpf(2) ** (-level)
    K = int(mp.mpf('4.6') * 2**level)
    S = mp.mpf(0); pih = mp.pi / 2
    for k in range(-K, K + 1):
        t = k * h
        sh = mp.sinh(t)
        u = mp.tanh(pih * sh)
        x = a + (b - a) * (1 + u) / 2
        dx = (b - a) / 2 * pih * mp.cosh(t) / (mp.cosh(pih * sh) ** 2)
        S += f(x) * dx
    return sgn * h * S

def J3val(a, b, c, p, q, level=5):
    if p == q: return mp.mpf(0)
    if a == p and b == a and c == a:
        return mp.log(abs(q - p))**3 / 6 * (1 if q > p else -1)
    lo, hi = min(p, q), max(p, q)
    inside = sorted(set(x for x in (b, c) if lo < x < hi and x != a))
    pts = [p] + inside + [q]
    tot = mp.mpf(0)
    if a == p:
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: J2val(b, c, p, t) / (t - a), u, v, level)
        return tot
    if a == q:
        fq = J2val(b, c, p, q)
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: (J2val(b, c, p, t) - fq) / (t - a) if t != a else mp.mpf(0), u, v, level)
        return tot + fq * mp.log(abs(q - p))
    if lo < a < hi:
        fa = J2val(b, c, p, a)
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: (J2val(b, c, p, t) - fa) / (t - a) if t != a else mp.mpf(0), u, v, level)
        return tot + fa * (mp.log(abs(q - a)) - mp.log(abs(a - p)))
    for u, v in zip(pts, pts[1:]):
        tot += ts(lambda t: J2val(b, c, p, t) / (t - a), u, v, level)
    return tot

_j3cache = {}
def J3(a, b, c, p, q):
    key = (a, b, c, p, q)
    if key not in _j3cache:
        _j3cache[key] = J3val(a, b, c, mp.mpf(p.numerator)/p.denominator if hasattr(p, 'numerator') else mp.mpf(p),
                                  mp.mpf(q.numerator)/q.denominator if hasattr(q, 'numerator') else mp.mpf(q))
    return _j3cache[key]

def delta_num(w, e, c):
    pts = (e,) + tuple(w) + (c,)
    tot = mp.mpf(0)
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = s05.lam2(gap, ci, ei)
                if not lam: continue
                tot += mp.mpf(lam.numerator) / lam.denominator * J3(ret[0], ret[1], ret[2], e, c)
            else:
                lam = s05.lam2(ret, c, e)
                if not lam: continue
                tot -= mp.mpf(lam.numerator) / lam.denominator * J3(gap[0], gap[1], gap[2], ei, ci)
    return tot

def delta_num_R4():
    tot = mp.mpf(0)
    for w, coef in s05.signed_words():
        tot += mp.mpf(coef.numerator) / coef.denominator * delta_num(w, 0, C0)
    return tot

# ---------- LLL fit ----------
l2, l3 = mp.log(2), mp.log(3)
pi2 = mp.pi**2
z3 = mp.zeta(3)
Li2a = mp.polylog(2, mp.mpf(1)/3)
L313, L31m = mp.polylog(3, mp.mpf(1)/3), mp.polylog(3, -mp.mpf(1)/2)
FITB = [('z3', z3), ('Li3(1/3)', L313), ('Li3(-1/2)', L31m),
        ('Li2*l2', Li2a*l2), ('Li2*l3', Li2a*l3),
        ('pi2*l2', pi2*l2), ('pi2*l3', pi2*l3),
        ('l2^3', l2**3), ('l2^2l3', l2**2*l3), ('l2*l3^2', l2*l3**2), ('l3^3', l3**3)]

def fit(tag, Dnum, precisions=(85, 105)):
    names = [nm for nm, _ in FITB]
    vals = [v for _, v in FITB]
    cands = {}
    for dps in precisions:
        mp.mp.dps = dps
        C = mp.mpf(10) ** (dps - 25)
        n = len(vals) + 1
        M = fmpz_mat(n, n + 1)
        for i, v in enumerate([Dnum] + vals):
            M[i, i] = fmpz(1)
            M[i, n] = fmpz(int(mp.nint(C * v)))
        R = M.lll(delta=0.99, eta=0.51)
        mp.mp.dps = DPS + 20
        row = next((r for r in [[R[r, cc] for cc in range(n + 1)] for r in range(n)] if r[0] != 0), None)
        coeffs = [int(cc) for cc in row[:-1]]
        g = reduce(gcd, [abs(cc) for cc in coeffs if cc]) or 1
        coeffs = [cc // g for cc in coeffs]
        if next((cc for cc in coeffs if cc), 1) < 0: coeffs = [-cc for cc in coeffs]
        cands[dps] = coeffs
        resid = sum(cc * v for cc, v in zip(coeffs, [Dnum] + vals))
        log(f"    {tag} fit dps {dps}: height {max(abs(cc) for cc in coeffs)}, resid@full {mp.nstr(abs(resid), 3)}")
    if cands[precisions[0]] != cands[precisions[1]]:
        log(f"    {tag}: NOT stable"); return None
    a = cands[precisions[1]]
    out = {}
    for cc, nm in zip(a[1:], names):
        if cc: out[nm] = Fr(-cc, a[0])
    log(f"    {tag}: " + ", ".join(f"{v}*{k}" for k, v in out.items()))
    return out

def delta_num_li32(x, y):
    return delta_num(s05.li32_word(x, y), 0, Fr(1))

if __name__ == '__main__':
    log(f"step_10 start ({time.strftime('%H:%M:%S')}), dps {DPS}")
    # controls
    c_li5 = delta_num((0, 0, 0, 0, 2), 0, C0)   # J = -Li5(1/4): classical
    log(f"control Li5(1/4): delta_num = {mp.nstr(c_li5, 5)} ({time.time()-t0:.0f}s); fit must give z3=Li3=0")
    fC = fit("Li5(1/4)-control", c_li5)
    if fC is not None:
        log(f"  control: z3={fC.get('z3', 0)}, Li3(1/3)={fC.get('Li3(1/3)', 0)}, Li3(-1/2)={fC.get('Li3(-1/2)', 0)}")
    D_R4 = delta_num_R4()
    log(f"D_num(R4^(5)(1/2)) computed ({time.time()-t0:.0f}s)")
    fR = fit("R4^(5)", D_R4)
    if fR:
        ok = fR.get('Li3(1/3)') == Fr(33, 2) and fR.get('Li3(-1/2)') == Fr(41)
        log(f"  Li3 coefficients match the EXACT symbol values (33/2, 41): {ok}")
        log(f"  a(R4^(5)) = {fR.get('z3')}")
    res = {}
    for (x, y) in [(Fr(1, 4), Fr(4, 3)), (Fr(8, 9), Fr(1))]:
        D_at = delta_num_li32(x, y)
        log(f"D_num(Li32({x},{y})) computed ({time.time()-t0:.0f}s)")
        res[(str(x), str(y))] = fit(f"Li32({x},{y})", D_at)
    if fR and all(res.values()):
        a_R = fR['z3']
        a_N5 = Fr(5, 4) * res[('1/4', '4/3')]['z3'] + Fr(1, 3) * res[('8/9', '1')]['z3']
        log(f"zeta3-cobracket: a(R4^(5)) = {a_R};  5/4*a(1/4,4/3) + 1/3*a(8/9,1) = {a_N5}")
        log(f"zeta3-part match (cobracket(R4^(5) - N5) = 0): {a_R == a_N5}")
    json.dump({'D_R4_fit': {k: str(v) for k, v in (fR or {}).items()},
               'control_fit': {k: str(v) for k, v in (fC or {}).items()},
               'atom_fits': {str(k): {k2: str(v2) for k2, v2 in v.items()} for k, v in res.items()}},
              open(G2B + 'zeta3_cobracket.json', 'w'), indent=1)
    log(f"step_10 done ({time.time()-t0:.0f}s)")
'''
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "step_10.py", "w").write(step10)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", "-n", "19", PY, "step_10.py"], cwd=G2B, capture_output=True, text=True, timeout=1800)
print("elapsed", time.time() - t0, "rc", r.returncode)
print(r.stdout)
print("STDERR:", r.stderr[-2500:])
