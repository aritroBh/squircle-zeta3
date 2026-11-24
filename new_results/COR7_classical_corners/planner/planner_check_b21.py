# Planner check of COR7 Task 2 (rational census) and Task 3 (c0 = sqrt2-1 closed form), independent of COR7's code.
# M3 from G1's proved B3 term lists (JSON) with the planner's own HPL series evaluator; weight split recomputed from the JSON.
import json, re, sys, math, mpmath as mp
from fractions import Fraction as Fr
from sympy import factorint, Matrix
D = int(sys.argv[1]) if len(sys.argv) > 1 else 330
mp.mp.dps = D + 20
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
C7 = "<project>/research/08_theory/COR7_classical_corners/"

# ---- Task 2: rational census by brute force (multiplicative rank of {2, c, 1-c, 1+c}) ----
def rank_rat(a, b):
    nums = [Fr(2), Fr(a, b), Fr(b - a, b), Fr(a + b, b)]
    primes = set()
    for q in nums:
        primes |= set(factorint(q.numerator)) | set(factorint(q.denominator))
    rows = []
    for q in nums:
        fn, fd = factorint(q.numerator), factorint(q.denominator)
        rows.append([fn.get(p, 0) - fd.get(p, 0) for p in sorted(primes)])
    return Matrix(rows).rank()
found = [(a, b) for b in range(2, 241) for a in range(1, b) if math.gcd(a, b) == 1 and rank_rat(a, b) <= 2]
print("rational c0 = a/b in (0,1), b <= 240, with rank <= 2:", found)

# ---- Task 3: M3 at c0 = sqrt2 - 1 ----
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
c = mp.sqrt(2) - 1
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
    part = mp.mpf(0)
    for kind, w, e, q in terms(G1 + fname):
        v = q * zval * term_val(kind, w, e)
        W[len(w) + zw] = W.get(len(w) + zw, 0) + v
        part += q * term_val(kind, w, e)
    if fname.endswith("(0,1,0).json"):
        Phi = mp.mpf(1)/2 - (1 - c**2)/(2*c) * mp.atanh(c)
        print("zeta(3) part + Phi(c) =", mp.nstr(part + Phi, 5), "(B3: zeta(3) part is -Phi)")
M3 = sum(W.values())
ref = mp.mpf(open(C7 + "b21_M3_hpl_1150.txt").read().split()[0])
print("weights present:", sorted(W))
print("M3(planner) =", mp.nstr(M3, 40))
print("|M3(planner) - COR7 b21_M3_hpl_1150| =", mp.nstr(abs(M3 - ref), 5))
# closed forms, transcribed from COR7 result.md Task 3
s2 = mp.sqrt(2); b = s2 - 1; L = mp.log(1 + s2); lam = mp.log(2); pi2 = mp.pi**2; z3 = mp.zeta(3)
Li = lambda n, x: mp.polylog(n, x)
W2 = -Li(2, b**2)/2 + pi2/12
W3 = ((-mp.mpf(5)/4 - 63*s2/32)*z3 + (-11 - 13*s2)*Li(3, b) + (5 + 4*s2)*Li(3, -b) + (14 + 16*s2)*Li(3, 1/s2)
      + (-4 - 4*s2)*L*Li(2, b) + (4 + 4*s2)*lam*Li(2, b) + (-mp.mpf(11)/12 - 49*s2/48)*pi2*L + (mp.mpf(29)/48 + 3*s2/4)*pi2*lam
      + (-mp.mpf(7)/24 - s2/3)*lam**3 + (-mp.mpf(7)/4 - 2*s2)*lam**2*L + (mp.mpf(11)/4 + 3*s2)*lam*L**2 + (-mp.mpf(1)/3 - s2/4)*L**3)
W4 = ((22 - 31*s2/2)*Li(4, b) + (24 - 27*s2/2)*Li(4, -b) + (-mp.mpf(5)/2 + 13*s2/8)*Li(4, -b**2) - s2*Li(4, 1/s2)
      + mp.mpf(3)/4*Li(4, mp.mpf(1)/2) + (-4 + 4*s2)*Li(4, 1 - 1/s2) + (-2 + 3*s2)*Li(4, 2 - s2)
      - 3*L*Li(3, b) + L*Li(3, -b) - 2*s2*L*Li(3, 1/s2) + (mp.mpf(53)/16 + 2*s2)*L*z3
      + 2*lam*Li(3, b) + (-2 + s2)*lam*Li(3, -b) + (-2 + 2*s2)*lam*Li(3, 1/s2) + (-mp.mpf(25)/16 - s2)*lam*z3
      + (2 + s2)*Li(2, b)**2 + (2 + s2/2)*Li(2, b)*L**2 + (-2 + s2/2)*Li(2, b)*lam*L + (mp.mpf(1)/4 - s2/4)*Li(2, b)*lam**2
      + (-mp.mpf(1)/2 - s2/8)*pi2*Li(2, b) + (mp.mpf(49)/960 - 29*s2/11520)*mp.pi**4
      + (-mp.mpf(5)/8 + s2/48)*pi2*L**2 + (mp.mpf(9)/16 - 11*s2/48)*pi2*lam*L + (-mp.mpf(3)/16 + 5*s2/32)*pi2*lam**2
      + (1 - s2/48)*L**4 + (-mp.mpf(5)/4 + s2/4)*L**3*lam + (mp.mpf(5)/8 - s2/4)*L**2*lam**2
      + (-mp.mpf(7)/24 + 13*s2/48)*L*lam**3 + (mp.mpf(13)/192 - 13*s2/384)*lam**4)
for w, cf in [(2, W2), (3, W3), (4, W4)]:
    print(f"W{w}: planner split = {mp.nstr(W[w], 25)}   |split - closed form| = {mp.nstr(abs(W[w] - cf), 5)}")
print("|M3 - (W2+W3+W4 closed forms)| =", mp.nstr(abs(M3 - W2 - W3 - W4), 5), " at dps", D)
