import json, re, sys, mpmath as mp
from fractions import Fraction as Fr
D = 300
mp.mp.dps = D + 20
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
WD = "<project>/research/08_theory/COR9_golden_corner_closed_forms/"

# ---- planner's independent HPL evaluator (verbatim from planner_check_b21.py) ----
def hpl(w, x, N):
    a = w[-1]; s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1): s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N): acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))
def terms(f):
    out = []
    for nm, v in json.load(open(f)):
        q = mp.mpf(Fr(v).numerator) / Fr(v).denominator
        if nm == '1': out.append(('1', (), 0, q)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        w = tuple(int(t) for t in m.group(3).split(','))
        out.append((m.group(2), w, int(m.group(1)), q))
    return out

for tag, c in [("b2sq", (mp.sqrt(2)-1)**2), ("gol1", 2/(1+mp.sqrt(5))), ("gol3", mp.sqrt(5)-2)]:
    N = int((D + 25) / -mp.log10(c)) + 30
    cache = {}
    def term_val(kind, w, e):
        if kind == '1': return mp.mpf(1)
        if (w, 1) not in cache: cache[(w, 1)] = hpl(w, c, N)
        f = c ** e * cache[(w, 1)]
        if kind == 'He':
            if (w, -1) not in cache: cache[(w, -1)] = hpl(w, -c, N)
            f += (-c) ** e * cache[(w, -1)]
        return f
    W = {}
    for fname, zw, zval in [("hpl_fit_k3_(0,0,0).json", 0, mp.mpf(1)), ("hpl_fit_k3_(1,0,0).json", 2, mp.zeta(2)), ("hpl_fit_k3_(0,1,0).json", 3, mp.zeta(3))]:
        for kind, w, e, q in terms(G1 + fname):
            W[len(w) + zw] = W.get(len(w) + zw, 0) + q * zval * term_val(kind, w, e)
    M3 = sum(W.values())
    ref = mp.mpf(open(WD+f"/{tag}_M3_hpl_1150.txt").read())
    print(f"{tag}: |planner-M3 - stored| = {mp.nstr(abs(M3-ref), 3)}  (weights {sorted(W)})")
