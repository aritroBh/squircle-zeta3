import sys, math
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import Li31
mp.mp.dps = 200
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
LHS = -Li31(mp.mpf(-1)/2, mp.mpf(-1))
RHS = 6*Li(4, mp.mpf(-2)) - 2*Li(4, mp.mpf(1)/2) - Li(3, mp.mpf(-1)/2)*l2 \
      - mp.mpf(3)/4*mp.zeta(3)*l2 + mp.mpf(5)/8*mp.pi**2*l2**2 + mp.mpf(137)/1440*mp.pi**4 + mp.mpf(1)/24*l2**4
print("LHS:", mp.nstr(LHS, 40))
print("RHS:", mp.nstr(RHS, 40))
print("|diff|:", mp.nstr(abs(LHS-RHS), 3))
# Let me ALSO compute Li4(-2) to 40 digits and double check via series
s = mp.nsum(lambda k: mp.mpf(-2)**k / k**4, [1, mp.inf])
print("Li4(-2):", mp.nstr(Li(4, mp.mpf(-2)), 40))
print("series:", mp.nstr(s, 40))
