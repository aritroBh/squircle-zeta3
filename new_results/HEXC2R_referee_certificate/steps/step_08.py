"""step_08.py -- referee TASK 4: recompute M3(2pi/3) from the COR5 closed form to >= 1100
digits with MY OWN Li2/Li3/Li4 (direct series for |z|<1; Jonquiere inversion for Li4(-2)),
compare against COR5's 1095-digit value and the hex_*.txt component values."""
import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)

DPS = 1200
mp.mp.dps = DPS

def my_polylog_series(s, z):
    """Li_s(z) for |z| < 1: absolutely convergent series, geometric tail bound."""
    z = mp.mpc(z)
    az = abs(z)
    assert az < 1
    eps = mp.mpf(10) ** (-DPS - 30)
    tot = mp.mpc(0); n = 1; zn = z
    while True:
        t = zn / mp.mpf(n) ** s
        tot += t
        # tail <= |t| * az/(1-az) once terms decrease geometrically (n^s monotone slow)
        if n > 100 and abs(t) * az / (1 - az) < eps:
            break
        n += 1; zn *= z
    return tot

def bernoulli4(y):
    """B_4(y) = y^4 - 2 y^3 + y^2 - 1/30."""
    return y**4 - 2*y**3 + y**2 - mp.mpf(1)/30

def my_li4_neg2():
    """Li_4(-2) by Jonquiere's inversion (my own):
    Li_4(-x) + Li_4(-1/x) = -(2 pi i)^4/4! B_4(1/2 + ln(x)/(2 pi i)), x = 2."""
    x = mp.mpf(2)
    y = mp.mpc(mp.mpf(1)/2) + mp.log(x) / (2 * mp.pi * mp.j)
    rhs = -(2 * mp.pi * mp.j) ** 4 / mp.factorial(4) * bernoulli4(y)
    return rhs - my_polylog_series(4, mp.mpf(-1)/2)

# cross-checks of my polylog vs mpmath
worst = 0
for s, z in [(2, "1/3"), (3, "1/3"), (3, "-1/2"), (4, "1/3"), (4, "-1/3"), (4, "2/3"),
             (4, "1/2"), (4, "-1/2"), (4, "1/9"), (2, "1/2"), (2, "-1/2")]:
    zf = mp.mpf(z)
    worst = max(worst, abs(my_polylog_series(s, zf) - mp.polylog(s, zf)))
log(f"my series Li_{{2,3,4}} vs mpmath at 11 points: worst {mp.nstr(worst, 3)}")
assert worst < mp.mpf("1e-1150")
li4m2_mine = my_li4_neg2()
li4m2_mp = mp.polylog(4, -2)
log(f"Li_4(-2): Jonquiere (mine) vs mpmath: |diff| = {mp.nstr(abs(li4m2_mine - li4m2_mp), 3)}")
assert abs(li4m2_mine - li4m2_mp) < mp.mpf("1e-1150")

pi, z2, z3 = mp.pi, mp.zeta(2), mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)
log(f"|zeta(2) - pi^2/6| = {mp.nstr(abs(z2 - pi**2/6), 3)} (the exact Euler certificate's identity)")
Li21 = my_polylog_series(2, mp.mpf(1)/3)
Li31 = my_polylog_series(3, mp.mpf(1)/3)
Li3m = my_polylog_series(3, mp.mpf(-1)/2)
Li41 = my_polylog_series(4, mp.mpf(1)/3)
Li4m3 = my_polylog_series(4, mp.mpf(-1)/3)
Li423 = my_polylog_series(4, mp.mpf(2)/3)
Li412 = my_polylog_series(4, mp.mpf(1)/2)
Li4m = my_polylog_series(4, mp.mpf(-1)/2)
Li419 = my_polylog_series(4, mp.mpf(1)/9)
Li4m2 = li4m2_mine

def rd(f):
    return mp.mpf(f.numerator) / f.denominator
from fractions import Fraction as Fr

W2 = L2**2/2 + Li21 + mp.log(mp.mpf(3)/2)**2/2
W3 = (Fr(-15, 4) and mp.mpf(-15)/4*z3) - 3*Li31 - 8*Li3m - pi**2*L2/12      + 4*Li21*L2 - mp.mpf(3)/2*Li21*L3 + L2**3/3 - mp.mpf(5)/2*L2**2*L3 + 2*L2*L3**2 - L3**3/4
F4 = (-21*Li4m2 - mp.mpf(65)/8*Li41 + mp.mpf(13)/2*Li423 + mp.mpf(137)/6*Li412 + mp.mpf(39)/64*Li419
      + mp.mpf(15)/2*Li31*L2 + mp.mpf(37)/2*Li3m*L2 - mp.mpf(5)/4*Li3m*L3
      + mp.mpf(119)/16*z3*L2 - mp.mpf(37)/32*z3*L3
      + Li21**2 - mp.mpf(23)/4*Li21*L2**2 + mp.mpf(7)/4*Li21*L2*L3 + Li21*L3**2 - pi**2*Li21/6
      - mp.mpf(455)/144*pi**2*L2**2 + mp.mpf(9)/8*pi**2*L2*L3 - mp.mpf(7)/32*pi**2*L3**2
      - mp.mpf(949)/2160*pi**4
      - mp.mpf(211)/144*L2**4 + mp.mpf(101)/24*L2**3*L3 - mp.mpf(15)/8*L2**2*L3**2
      - mp.mpf(35)/24*L2*L3**3 + mp.mpf(113)/192*L3**4)
# Psi2(1/2) closed: from the fit: (3/8)(H01+H0-1) + (1/2)H_{1,-1}... use COR5-consistent value:
# my own HPL evaluation of Psi2 (step in referee notes): recompute here with my HPL series
sys.path.insert(0, "../COR4_M3alpha_classical")
import json, re
import reflib as R
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
psi2 = mp.mpf(0)
for key, val in json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json")):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
    if not m: continue
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    if len(w) != 2: continue
    psi2 += mp.mpf(str(val)) * mp.mpf(1)/2**e * R.hpl_value(w, mp.mpf(1)/2)
W4 = F4 + z2 * psi2 + mp.mpf(3)/8 * z3 * L3
M3 = W4 + W3 + W2

C5 = "../COR5_corner_classicality/"
def ref(name):
    return mp.mpf(open(C5 + name).read().strip())
checks = [
    ("W2", W2, "hex_W2.txt"), ("W3", W3, "hex_W3.txt"), ("W4", W4, "hex_W4.txt"),
    ("F4", F4, "hex_F4.txt"), ("Psi2", psi2, "hex_Psi2.txt"),
]
for name, v, f in checks:
    d = abs(v - ref(f))
    log(f"{name}: |closed_form - {f}| = {mp.nstr(d, 3)}")
    assert d < mp.mpf("1e-1080")
m3ref = open(C5 + "M3_2pi3_1100.txt").read().strip()
diff = abs(M3 - mp.mpf(m3ref))
log(f"M3(2pi/3) closed form = {mp.nstr(M3, 60)}")
log(f"|M3 closed form - COR5 1095-digit value| = {mp.nstr(diff, 3)}")
ndig = int(-mp.log10(diff)) if diff > 0 else 9999
log(f"agreement: {ndig} digits (reference has ~1095)")
assert diff < mp.mpf("1e-1090")
log("STEP 08 PASS")
with open("progress.log", "a") as f:
    f.write("step_08: M3(2pi/3) closed form recomputed at 1200 dps, agrees %d digits with COR5 1095-digit value\n" % ndig)
