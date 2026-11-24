import sys, math
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series, Li31
mp.mp.dps = 1450
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
N = int((1450 + 15) / abs(math.log10(0.5))) + 10
c = mp.mpf(1)/2
Hv = H_series((0,0,-1,1), c, N)
RHS = 6*Li(4, mp.mpf(-2)) - 2*Li(4, c) - Li(3, -c)*l2 - mp.mpf(3)/4*mp.zeta(3)*l2 \
      + mp.mpf(5)/8*mp.pi**2*l2**2 + mp.mpf(137)/1440*mp.pi**4 + mp.mpf(1)/24*l2**4
print("dps 1450: |RHS - Hv| =", mp.nstr(abs(RHS - Hv), 3))
# and verify Hv = -Li31(-1/2,-1)
print("|Hv + Li31(-1/2,-1)| =", mp.nstr(abs(Hv + Li31(-c, -mp.mpf(1))), 3))
# ALSO check: Li4(-2) value vs a series: Li4(-2) = sum (-2)^k/k^4
s = mp.nsum(lambda k: mp.mpf(-2)**k / k**4, [1, mp.inf])
print("Li4(-2) series check:", mp.nstr(abs(s - Li(4, mp.mpf(-2))), 3))
print("Hv:", mp.nstr(Hv, 25))
