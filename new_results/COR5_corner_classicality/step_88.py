import sys, math
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series, Li31
mp.mp.dps = 1120
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
N = int((1120 + 15) / abs(math.log10(0.5))) + 10
c = mp.mpf(1)/2
Hv = H_series((0,0,-1,1), c, N)
print("H(0,0,-1,1)(1/2) =", mp.nstr(Hv, 30))
print("vs -Li31(-1/2,-1):", mp.nstr(abs(Hv - (-Li31(-c, -mp.mpf(1)))), 3))
RHS = 6*Li(4, mp.mpf(-2)) - 2*Li(4, c) - Li(3, -c)*l2 - mp.mpf(3)/4*mp.zeta(3)*l2 \
      + mp.mpf(5)/8*mp.pi**2*l2**2 + mp.mpf(137)/1440*mp.pi**4 + mp.mpf(1)/24*l2**4
print("|RHS - Hv| =", mp.nstr(abs(RHS - Hv), 3))
# Also: symbol-predicted combination Hv - 3*Li4(-2) reduce? check numeric value
# what IS Hv's class then? Let me check the pipeline symbol quotient coords of alternative
# maybe the symbol of H(0,0,-1,1) should be computed as I(0;1,-1,0,0): FIRST letter omega_1...
# try the REVERSED convention: S = g(a1)...g(an) (no reversal):
import sys
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import W4Quotient
Q = W4Quotient(2)
e1 = (Fr(1), Fr(0))
word_alt = ((Fr(-1),Fr(0)), (Fr(-1),Fr(0)), (Fr(-1),Fr(1)), (Fr(1),Fr(0)))  # g(0)g(0)g(-1)g(1) unreversed
print("alt-convention coords:", Q.coords({word_alt: Fr(1)}))
