"""verify step: from OWN values (verify_vals.py) fit the polynomial in k of W_{d+1}(k), k >= d+1 (degree d-1),
check the remaining k, and compare the top coefficient with the claimed law  g_m = eps_m A_{m+1}/(2 m!) [y^{m+2}] E(y)
(m = d-1) and, where claimed, with the closed forms for c_4 (weight 6), c_2^(4), c_1^(4) (weight 4), c_3^(5) (weight 5).
usage: verify_fit.py OWNDIR DPS D
Prints mp.mp.dps next to every residual.  mp.mp.dps is set AFTER all imports.
"""
import sys, os, json, glob
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp

OWN = sys.argv[1]; DPS = int(sys.argv[2]); D = int(sys.argv[3])
mp.mp.dps = DPS + 10          # AFTER imports
print(f"mp.mp.dps = {mp.mp.dps}  (data target {DPS} digits)", flush=True)

# ------------------------------------------------------------ data
vals = {}
for fn in glob.glob(os.path.join(OWN, f"own_k*_d{D}_{DPS}.json")):
    r = json.load(open(fn))
    vals[r["k"]] = Fr(r["value"])
ks = sorted(vals)
print(f"d={D} (weight {D+1}); own data at k =", ks)
deg = D - 1
k0 = D + 1                       # first k of the regime (k >= weight)
reg = [k for k in ks if k >= k0]
assert len(reg) >= deg + 1, "not enough points"

def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)

def poly_through(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co

fitk = reg[:deg + 1]
co = poly_through([(Fr(k), vals[k]) for k in fitk])
print("polynomial of degree", deg, "through k =", fitk)
for j, c in enumerate(co):
    print(f"   c_{j} = {mp.nstr(tomp(c), 60)}")
for k in ks:
    pk = sum(co[j] * Fr(k) ** j for j in range(len(co)))
    tag = "fitted" if k in fitk else ("check" if k >= k0 else "BELOW regime")
    print(f"   k={k:2d} [{tag}]: W - p(k) = {mp.nstr(tomp(vals[k] - pk), 5)}  (mp.mp.dps={mp.mp.dps})")
# finite differences
def diffs(seq, n):
    cur = list(seq)
    for _ in range(n):
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur
ys = [vals[k] for k in reg]
for n in (deg, deg + 1):
    if len(ys) > n:
        dn = diffs(ys, n)
        print(f"   regime k={reg[0]}..{reg[-1]}: max |Delta^{n}| = {mp.nstr(max(abs(tomp(x)) for x in dn), 5)}")

# ------------------------------------------------------------ the claimed law for the leading coefficient
# E(y) = exp(sum_n s_n y^n / n),  s_n = eta(n) (n odd), zeta(n) (n even); eta(1) = log 2
NN = 12
s = [None] + [(mp.log(2) if n == 1 else (mp.altzeta(n) if n % 2 == 1 else mp.zeta(n))) for n in range(1, NN + 1)]
L1 = [mp.mpf(0)] + [s[n] / n for n in range(1, NN + 1)]      # series of log E
Ec = [mp.mpf(1)] + [mp.mpf(0)] * NN                          # exp via  n E_n = sum_{j} j L_j E_{n-j}
for n in range(1, NN + 1):
    Ec[n] = sum(j * L1[j] * Ec[n - j] for j in range(1, n + 1)) / n
# zigzag (Euler up/down) numbers A_n: 1,1,1,2,5,16,61,272,1385,...
def zigzag(n):
    A = [1]
    # Seidel-Entringer / boustrophedon: A_n = number of alternating permutations
    row = [1]
    out = [1]
    for i in range(1, n + 1):
        new = [0]
        for x in reversed(row):
            new.append(new[-1] + x)
        row = new
        out.append(row[-1])
    return out
Az = zigzag(10)
assert Az[:9] == [1, 1, 1, 2, 5, 16, 61, 272, 1385], Az
eps = lambda m: [1, -1, -1, 1][m % 4]
m = D - 1
g = eps(m) * Az[m + 1] * Ec[m + 2] / (2 * mp.factorial(m))
print(f"law g_{m} = eps_m A_(m+1)/(2 m!) [y^{m+2}]E = {mp.nstr(g, 60)}")
print(f"   fitted leading coefficient c_{deg} - law = {mp.nstr(tomp(co[deg]) - g, 5)}   (mp.mp.dps={mp.mp.dps})")
# independent evaluation of [y^(m+2)] E(y) from the Gamma-function form by numerical Taylor expansion
if mp.mp.dps <= 160:
    Ef = lambda y: mp.mpf(2) ** y * mp.gamma(1 - y) * mp.gamma(1 + y / 2) / mp.gamma(1 - y / 2)
    tay = mp.taylor(Ef, 0, m + 2)
    print(f"   [y^{m+2}]E from Gamma form (mp.taylor) - series = {mp.nstr(tay[m+2] - Ec[m+2], 5)}")

# ------------------------------------------------------------ closed forms
lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi = mp.pi; pi2 = pi ** 2
r2 = mp.sqrt(2); b = r2 - 1; Lp = mp.log(1 + r2)
if D == 5:
    c4 = (168 * lam ** 6 + 420 * pi2 * lam ** 4 + 378 * pi2 ** 2 * lam ** 2 + 61 * pi2 ** 3 + 5040 * z3 * lam ** 3
          + 2520 * pi2 * z3 * lam + 3780 * z3 ** 2 + 22680 * z5 * lam) / 362880
    print(f"c_4 closed form (PSLQ relation) = {mp.nstr(c4, 60)}")
    print(f"   fitted c_4 - closed form = {mp.nstr(tomp(co[4]) - c4, 5)}   (mp.mp.dps={mp.mp.dps})")
    print(f"   law g_4 - closed form    = {mp.nstr(g - c4, 5)}")
if D == 3:
    c2 = -(120 * z3 * lam + 3 * pi ** 4 + 20 * pi2 * lam ** 2 + 20 * lam ** 4) / 960
    li3 = lambda x: mp.polylog(3, x)
    c1 = (11 * pi ** 4 - 210 * Lp * z3 - 100 * pi2 * Lp ** 2 + 160 * Lp ** 4 + 360 * lam * z3 + 100 * pi2 * lam * Lp
          + 80 * lam * Lp ** 3 + 55 * pi2 * lam ** 2 - 120 * lam ** 2 * Lp ** 2 - 40 * lam ** 3 * Lp + 40 * lam ** 4
          - 480 * lam * (li3(b) + li3(-b)) - 960 * Lp * (li3(b) - li3(-b)) + 1920 * Lp * li3(1 / r2)) / 480
    print(f"c_2^(4) closed form = {mp.nstr(c2, 50)};  fitted - closed = {mp.nstr(tomp(co[2]) - c2, 5)}  (mp.mp.dps={mp.mp.dps})")
    print(f"c_1^(4) closed form = {mp.nstr(c1, 50)};  fitted - closed = {mp.nstr(tomp(co[1]) - c1, 5)}  (mp.mp.dps={mp.mp.dps})")
if D == 4:
    c3 = (270 * z5 + 30 * pi2 * z3 + 180 * z3 * lam ** 2 + 9 * pi ** 4 * lam + 20 * pi2 * lam ** 3 + 12 * lam ** 5) / 3456
    print(f"c_3^(5) closed form = {mp.nstr(c3, 50)};  fitted - closed = {mp.nstr(tomp(co[3]) - c3, 5)}  (mp.mp.dps={mp.mp.dps})")
if D == 5 and 5 in vals:
    pass
