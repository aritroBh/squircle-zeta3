"""Independent check of M3 = W4 + W3' + T: W4 by direct quadrature of its integral (eq. W4 of the note)."""
import os, sys, pickle
import mpmath as mp
import polylog_exact  # noqa: F401  (import before setting precision: dilog_exact sets dps=60 on import)

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 120
mp.mp.dps = DPS + 20
r2 = mp.sqrt(2); a = 1 + r2; b = r2 - 1
l1 = lambda v: mp.log((v - b) * (v + a) / (4 * v))
l2 = lambda v: mp.log((a - v) * (v + b) / (4 * v))
f = lambda v: l1(v) ** 3 * ((8 - 7 * r2) / (48 * v) - 1 / (3 * (v + a))) + r2 / 16 * l1(v) ** 2 * l2(v) / v
W4q = mp.quad(f, [b, mp.mpf(1) / 2, 1, 2, a])
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs")
t = pickle.load(open(os.path.join(OUT, "m3_tails_exact.pkl"), "rb"))
M3 = mp.mpf(open(os.path.join(OUT, "M3_methodA_300.txt")).read().strip())
W4s = M3 - t["W3p"].num() - t["T"].num()
with mp.workdps(DPS):
    print("W4 by quadrature       =", mp.nstr(W4q, 60))
    print("W4 = M3 - W3' - T      =", mp.nstr(W4s, 60))
    print("difference             =", mp.nstr(W4q - W4s, 3))
