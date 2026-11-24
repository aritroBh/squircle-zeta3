import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 1120
M3 = mp.mpf(open(WD+"b2sq_M3_hpl_1150.txt").read())
W4v = mp.mpf(open(WD+"b2sq_W4_hpl_1150.txt").read())
W3v = mp.mpf(open(WD+"b2sq_W3_hpl_1150.txt").read())
W2v = mp.mpf(open(WD+"b2sq_W2_hpl_1150.txt").read())
b = mp.sqrt(2)-1; rt2 = mp.sqrt(2); c = b*b
Llog = mp.log(1+rt2); lam = mp.log(2)
Li = lambda s, x: mp.polylog(s, x)
z3, pi2, pi4 = mp.zeta(3), mp.pi**2, mp.pi**4
Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2)
Li2b = Li(2, b)
r4 = json.load(open(WD+"b2sq_W4_relation.json")); a0, cf = r4["a0"], r4["coeffs"]
B4 = [Li(4, b), Li(4, -b), Li(4, -b*b), Li(4, 1/rt2), Li(4, mp.mpf(1)/2), Li(4, 1-1/rt2), Li(4, 2-rt2),
      Llog*Li3b, Llog*Li3mb, Llog*Li3s, Llog*z3,
      lam*Li3b, lam*Li3mb, lam*Li3s, lam*z3,
      Li2b**2, Li2b*Llog**2, Li2b*lam*Llog, Li2b*lam**2, pi2*Li2b,
      pi4, pi2*Llog**2, pi2*lam*Llog, pi2*lam**2,
      Llog**4, Llog**3*lam, Llog**2*lam**2, Llog*lam**3, lam**4]
W4cf = -sum(mp.mpf(cf[i])*B4[i] + mp.mpf(cf[i+29])*rt2*B4[i] for i in range(29))/a0
r3 = json.load(open(WD+"b2sq_W3_relation.json")); a03, cf3 = r3["a0"], r3["coeffs"]
B3 = [z3, Li3b, Li3mb, Li3s, Llog*Li2b, lam*Li2b, pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
W3cf = -sum(mp.mpf(cf3[i])*B3[i] + mp.mpf(cf3[i+12])*rt2*B3[i] for i in range(12))/a03
W2cf = -Li(2, c) - Li(2, -c) + pi2/12
print("|W4 - W4cf| =", mp.nstr(abs(W4v - W4cf), 3))
print("|W3 - W3cf| =", mp.nstr(abs(W3v - W3cf), 3))
print("|W2 - W2cf| =", mp.nstr(abs(W2v - W2cf), 3))
print("|M3 - (W4cf+W3cf+W2cf)| =", mp.nstr(abs(M3 - (W4cf + W3cf + W2cf)), 3))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_30: b2sq assembly; Task 1 complete.\n")
