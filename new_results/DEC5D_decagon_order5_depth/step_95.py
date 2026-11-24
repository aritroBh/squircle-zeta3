
import pickle
import mpmath as mp
import flint
from flint import fmpz_mat
mp.mp.dps = 420
# L(chi_5, 5) via Hurwitz zeta
chi = {1:1, 2:-1, 3:-1, 4:1}
Lchi5_5 = sum(mp.mpf(c)*mp.zeta(5, mp.mpf(a)/5) for a, c in chi.items())/5**5
print("L(chi_5,5) =", mp.nstr(Lchi5_5, 50))
# sanity vs direct series
s = mp.nsum(lambda n: [-1,1,-1,-1,1][n%5]/mp.mpf(n)**5, [1, mp.inf])
print("direct series:", mp.nstr(s, 50), " diff:", mp.nstr(abs(Lchi5_5-s),3))
