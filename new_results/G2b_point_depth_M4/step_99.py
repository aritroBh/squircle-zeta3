
step10 = r'''"""step_10.py -- G2b: the numerical (3,2) cobracket at the hexagon, including the zeta(3)-part
invisible to the symbol.

delta_word (M4c's formula, ported and exact-checked at the symbol level in step_05) with:
  - exact Lambda^2 coordinates lam (step_05 machinery),
  - retained 3-letter from-0 integrals by HPL series (exact to working precision),
  - k=3 gap integrals J3(a,b,c; p->q) by closed-form J2 + one-level quadrature:
      J2(b,c;p,t) closed form (validated vs series 1e-102 and vs quadrature 1e-102),
      PV subtraction at interior first letters, tangential regularization at endpoints.
D_num(P) is then fitted over {zeta3, Li3(1/3), Li3(-1/2), weight-3 products} by LLL:
the Li3 coefficients must (and do) match the EXACT symbol values (validation of the whole
numerical chain); the zeta(3) coefficient a_P is the new information.

Controls: delta_num(Li5(z)) ~ 0 and delta_num(shuffle products) ~ 0 (to quad precision).
Then: a(R4^(5)) vs 5/4*a(Li32(1/4,4/3)) + 1/3*a(Li32(8/9,1))  -- the zeta(3)-cobracket match.
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

DPS = 260
mp.mp.dps = DPS + 30
t0 = time.time()
C0 = Fr(1, 2)

def log(msg):
    print(msg, flush=True)
    with open(G2B + "progress.log", "a") as f:
        f.write(msg + "\n")

# ---------- closed-form J2 and one-level J3 (validated) ----------
def Gfun(u, b, c):
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + mp.polylog(2, (u - c) / (b - c)).real

def J2val(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b:
            return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

def J3val(a, b, c, p, q):
    """regularized/PV J((a,b,c); p->q) = int_p^q J2(b,c;p,t)/(t-a) dt."""
    if p == q: return mp.mpf(0)
    lo, hi = min(p, q), max(p, q)
    inside = sorted(set(x for x in (b, c) if lo < x < hi and x != a))
    pts = [p] + inside + [q]
    if a == p:
        f = lambda t: J2val(b, c, p, t) / (t - a)
        return mp.quad(f, pts)
    if a == q:
        fq = J2val(b, c, p, q)
        g = lambda t: (J2val(b, c, p, t) - fq) / (t - a) if t != a else mp.mpf(0)
        return mp.quad(g, pts) + fq * mp.log(abs(q - p))
    if lo < a < hi:
        fa = J2val(b, c, p, a)
        g = lambda t: (J2val(b, c, p, t) - fa) / (t - a) if t != a else mp.mpf(0)
        pv = fa * (mp.log(abs(q - a)) - mp.log(abs(a - p)))
        return mp.quad(g, pts) + pv
    f = lambda t: J2val(b, c, p, t) / (t - a)
    return mp.quad(f, pts)

# ---------- from-0 3-letter integrals via HPL series ----------
def hpl_series_val(w, x, N):
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))

NSER = int((DPS + 30) / mp.log10(2)) + 80
from functools import lru_cache
@lru_cache(maxsize=None)
def Jval0(w, x):   # J(w; 0->x), letters {0,+-1}, x in {+-1, 1/2}: J = (-1)^#1 H_w(x)
    return ((-1) ** sum(1 for a in w if a == 1)) * hpl_series_val(w, mp.mpf(x), NSER)

def delta_num(w, e, c):
    """numerical cobracket value (coefficient of Li2(1/3) in the Lambda^2 sense)."""
    pts = (e,) + tuple(w) + (c,)
    tot = mp.mpf(0)
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = s05.lam2(gap, ci, ei)
                if not lam: continue
                tot += mp.mpf(lam.numerator) / lam.denominator * Jval0(ret, c)
            else:
                lam = s05.lam2(ret, c, e)
                if not lam: continue
                a3, b3, c3 = gap
                tot -= mp.mpf(lam.numerator) / lam.denominator * J3val(a3, b3, c3, ei, ci)
    return tot

def delta_num_R4():
    tot = mp.mpf(0)
    for w, coef in s05.signed_words():
        tot += mp.mpf(coef.numerator) / coef.denominator * delta_num(w, 0, C0)
    return tot

# ---------- LLL fit over {z3, Li3(1/3), Li3(-1/2), products} ----------
l2, l3 = mp.log(2), mp.log(3)
pi2 = mp.pi**2
z3 = mp.zeta(3)
Li2a = mp.polylog(2, mp.mpf(1)/3)
L313, L31m = mp.polylog(3, mp.mpf(1)/3), mp.polylog(3, -mp.mpf(1)/2)
FITB = [('z3', z3), ('Li3(1/3)', L313), ('Li3(-1/2)', L31m),
        ('Li2*l2', Li2a*l2), ('Li2*l3', Li2a*l3),
        ('pi2*l2', pi2*l2), ('pi2*l3', pi2*l3),
        ('l2^3', l2**3), ('l2^2l3', l2**2*l3), ('l2*l3^2', l2*l3**2), ('l3^3', l3**3)]

def fit(tag, Dnum, precisions=(190, 240)):
    names = [nm for nm, _ in FITB]
    vals = [v for _, v in FITB]
    cands = {}
    for dps in precisions:
        mp.mp.dps = dps
        C = mp.mpf(10) ** (dps - 30)
        n = len(vals) + 1
        M = fmpz_mat(n, n + 1)
        for i, v in enumerate([Dnum] + vals):
            M[i, i] = fmpz(1)
            M[i, n] = fmpz(int(mp.nint(C * v)))
        R = M.lll(delta=0.99, eta=0.51)
        mp.mp.dps = DPS + 30
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
    w = s05.li32_word(x, y)
    return delta_num(w, 0, Fr(1))

if __name__ == '__main__':
    log(f"step_10 start ({time.strftime('%H:%M:%S')}), dps {DPS}")
    # controls: Li5 and shuffle products must have zero cobracket value
    zt = Fr(1, 3)
    w_li5 = (0, 0, 0, 0, 3)      # J = -Li5(1/3)
    c_li5 = delta_num(w_li5, 0, C0)
    log(f"control Li5(1/3): delta_num = {mp.nstr(c_li5, 5)} (must be ~0)")
    # shuffle control: J((0,1);0->1/2)*J((0,0,1);0->1/2) as shuffle sum of words
    from cor6_sym import shuffles as _sh
    dsh = mp.mpf(0)
    for w in _sh((0, 1), (0, 0, 1)):
        dsh += delta_num(w, 0, C0)
    log(f"control shuffle (0,1)x(0,0,1): delta_num = {mp.nstr(dsh, 5)} (must be ~0)")
    # the target
    D_R4 = delta_num_R4()
    log(f"D_num(R4^(5)(1/2)) computed ({time.time()-t0:.0f}s)")
    fR = fit("R4^(5)", D_R4)
    if fR:
        ok = fR.get('Li3(1/3)') == Fr(33, 2) and fR.get('Li3(-1/2)') == Fr(41)
        log(f"  Li3 coefficients match the EXACT symbol values (33/2, 41): {ok}")
        log(f"  a(R4^(5)) = {fR.get('z3')}")
    # atoms
    res = {}
    for (x, y) in [(Fr(1, 4), Fr(4, 3)), (Fr(8, 9), Fr(1))]:
        D_at = delta_num_li32(x, y)
        fa = fit(f"Li32({x},{y})", D_at)
        res[(str(x), str(y))] = fa
    if fR and all(res.values()):
        a_R = fR['z3']
        a_N5 = Fr(5, 4) * res[('1/4', '4/3')]['z3'] + Fr(1, 3) * res[('8/9', '1')]['z3']
        log(f"zeta3-cobracket: a(R4^(5)) = {a_R}, 5/4*a(Li32(1/4,4/3)) + 1/3*a(Li32(8/9,1)) = {a_N5}")
        log(f"zeta3-part match: {a_R == a_N5}")
    json.dump({'D_R4_fit': {k: str(v) for k, v in (fR or {}).items()},
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
