import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 1120
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2
Li = lambda s, x: mp.polylog(s, x)
for tag, c in [("gol1", 1/phi), ("gol3", r5-2)]:
    W2v = mp.mpf(open(WD+f"{tag}_W2_hpl_1150.txt").read())
    W2a = -Li(2, c) - Li(2, -c) + pi2/12          # direct from HPL weight-2 terms
    W2b = -Li(2, c*c)/2 + pi2/12                  # duplication
    print(f"{tag}: |W2 - (-Li2(c)-Li2(-c)+pi^2/12)| = {mp.nstr(abs(W2v-W2a),3)}   |W2 - (-Li2(c^2)/2+pi^2/12)| = {mp.nstr(abs(W2v-W2b),3)}")
# gol1 elementary form: pi^2/20 + L^2/2
W2v = mp.mpf(open(WD+"gol1_W2_hpl_1150.txt").read())
print("gol1: |W2 - (pi^2/20 + L^2/2)| =", mp.nstr(abs(W2v - (pi2/20 + L**2/2)), 3))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_47: golden W2 closed forms.\n")
