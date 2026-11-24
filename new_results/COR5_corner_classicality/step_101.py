import sys, math
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import Li31
mp.mp.dps = 2100
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
LHS = -Li31(mp.mpf(-1)/2, mp.mpf(-1))
RHS = 6*Li(4, mp.mpf(-2)) - 2*Li(4, mp.mpf(1)/2) - Li(3, mp.mpf(-1)/2)*l2 \
      - mp.mpf(3)/4*mp.zeta(3)*l2 + mp.mpf(5)/8*mp.pi**2*l2**2 + mp.mpf(137)/1440*mp.pi**4 + mp.mpf(1)/24*l2**4
d = abs(LHS - RHS)
print("dps 2100: |LHS-RHS| =", mp.nstr(d, 3))
print("LHS:", mp.nstr(LHS, 50))
print("RHS:", mp.nstr(RHS, 50))
