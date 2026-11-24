
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
verify_g2 = r'''"""verify_g2.py -- < 5 min re-check of G2's B4 at two angles (the one G2 did not deliver).

Adapted from planner_check_b4.py (independent mpmath HPL series evaluator, Remiddi-Vermaseren):
  angle pi/2 (square, c = 1/sqrt(2)): |B4 - M4_800| (M4a full-integrand quadrature, 800 digits);
  angle 2pi/3 (hexagon, c = 1/2):     |B4 - G2's values4_2pi_3.txt quadrature| (~92-digit limit).
Runs at dps 110 (planner's default).
"""
import json, re, sys, mpmath as mp
from fractions import Fraction as Fr
D = 110
mp.mp.dps = D + 20
G2 = "<project>/research/08_theory/G2_M4alpha/"
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"

def hpl(w, x, N):
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

R4 = terms(G2 + "hpl_fit_k4_(0,0,0).json")
P2 = terms(G2 + "hpl_fit_k4_(1,0,0).json")
PS = terms(G1 + "hpl_fit_k3_(1,0,0).json")

def B4(c):
    N = int((D + 25) / -mp.log10(abs(c))) + 30
    cache = {}
    at = mp.atanh(c)
    F1 = (7*c/16 - mp.mpf(5)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(9)/16 - mp.mpf(3)/(16*c**2)
    F2 = (-c/16 - mp.mpf(1)/(8*c) + mp.mpf(3)/(16*c**3))*at + mp.mpf(1)/16 - mp.mpf(3)/(16*c**2)
    return (ev(R4, c, N, cache) + mp.zeta(2)*ev(P2, c, N, cache) - mp.zeta(3)*ev(PS, c, N, cache)
            + mp.zeta(4)*F1 + mp.zeta(2)**2*F2)

ok = True
M4sq = mp.mpf(open("<project>/research/08_theory/M4a_decomposition/M4_800.txt").read().split()[0])
v = B4(1/mp.sqrt(2))
r1 = abs(v - M4sq)
print(f"square corner c=1/sqrt2: |B4 - M4_800| = {mp.nstr(r1, 5)}")
ok &= r1 < mp.mpf(10) ** (-D + 15)
qline = [l for l in open(G2 + "values4_2pi_3.txt") if l.startswith("M4 =")][0]
v = B4(mp.cos(mp.pi/3))
r2 = abs(v - mp.mpf(qline.split("=")[1].strip()))
print(f"hexagon corner c=1/2 (2pi/3): |B4 - quadrature| = {mp.nstr(r2, 5)}  (quadrature limit ~1e-92)")
ok &= r2 < mp.mpf(10) ** (-90)
print(f"verify_g2: {'ALL PASS' if ok else 'FAIL'} (dps {D})")
sys.exit(0 if ok else 1)
'''
open(G2B + "verify_g2.py", "w").write(verify_g2)
import subprocess, time
PY = "<project>/.venv/bin/python"
for f in ["verify_g2b.py", "verify_g2.py"]:
    t0 = time.time()
    r = subprocess.run(["nice", "-n", "19", PY, f], cwd=G2B, capture_output=True, text=True, timeout=300)
    print(f"===== {f} ({time.time()-t0:.0f}s, rc {r.returncode}) =====")
    print(r.stdout)
    print("STDERR:", r.stderr[-1200:])
