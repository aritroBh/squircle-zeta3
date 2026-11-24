# Planner check of G2's B4 (M4 at every angle), with an independent mpmath HPL series evaluator.
# Uses only G2's/G1's term lists (JSON) and the F1, F2 formulas printed in G2's result.md.
import json, re, sys, mpmath as mp
from fractions import Fraction as Fr
D = int(sys.argv[1]) if len(sys.argv) > 1 else 110
mp.mp.dps = D + 20
G2 = "<project>/research/08_theory/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"

def hpl(w, x, N):
    """H_w(x), letters 0,+-1, w not ending in 0: power series to N terms (Remiddi-Vermaseren)."""
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

def terms(f):
    out = []
    for nm, v in json.load(open(f)):
        q = mp.mpf(Fr(v).numerator) / Fr(v).denominator
        if nm == '1': out.append(('1', None, 0, q)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        if not m: raise ValueError(nm)
        w = tuple(int(t) for t in m.group(3).split(','))
        out.append((m.group(2), w, int(m.group(1)), q))
    return out

def ev(ts, c, N, cache):
    tot = mp.mpf(0)
    for kind, w, e, q in ts:
        if kind == '1': tot += q; continue
        if (w, c) not in cache: cache[(w, c)] = hpl(w, c, N)
        f = c ** e * cache[(w, c)]
        if kind == 'He':
            if (w, -c) not in cache: cache[(w, -c)] = hpl(w, -c, N)
            f += (-c) ** e * cache[(w, -c)]
        tot += q * f
    return tot

R4 = terms(G2 + "hpl_fit_k4_(0,0,0).json"); P2 = terms(G2 + "hpl_fit_k4_(1,0,0).json"); PS = terms(G1 + "hpl_fit_k3_(1,0,0).json")
def B4(c):
    N = int((D + 25) / -mp.log10(abs(c))) + 30
    cache = {}
    at = mp.atanh(c)
    F1 = (7*c/16 - mp.mpf(5)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(9)/16 - mp.mpf(3)/(16*c**2)
    F2 = (-c/16 - mp.mpf(1)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(1)/16 - mp.mpf(3)/(16*c**2)
    return (ev(R4, c, N, cache) + mp.zeta(2)*ev(P2, c, N, cache) - mp.zeta(3)*ev(PS, c, N, cache)
            + mp.zeta(4)*F1 + mp.zeta(2)**2*F2)

M4sq = mp.mpf(open("<project>/research/08_theory/M4a_decomposition/M4_800.txt").read().split()[0])
v = B4(1/mp.sqrt(2)); print("square corner c=1/sqrt2: |B4 - M4_800| =", mp.nstr(abs(v - M4sq), 5))
def qval(f):
    for line in open(G2 + f):
        if line.startswith("M4 ="): return mp.mpf(line.split("=")[1].strip())
for f, c in [("values4_2pi_3.txt", mp.cos(mp.pi/3)), ("values4_pi_3.txt", mp.cos(mp.pi/6))]:
    v = B4(c); print(f"{f}: c={mp.nstr(c,6)}: |B4 - quadrature| =", mp.nstr(abs(v - qval(f)), 5))
print("dps", D)
