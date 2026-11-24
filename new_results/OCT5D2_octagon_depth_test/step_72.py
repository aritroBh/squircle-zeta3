import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 250
from fractions import Fraction as F
from flint import fmpz_mat

r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2); b = r2-1
def Li32(x, y, N=1400):
    tot = mp.mpf(0); xp = x
    for n in range(2, N+1):
        if n > 2: xp *= x
        inner = mp.mpf(0); yp = y
        for m in range(1, n):
            if m > 1: yp *= y
            inner += yp/(m*m)
        tot += xp/(n**3)*inner
    return tot
z5 = mp.zeta(5); z3 = mp.zeta(3); z2 = mp.pi**2/6
t0=time.time()
z32 = Li32(mp.mpf(1), mp.mpf(1), 900)   # zeta(3,2), slow-ish convergence: use integral instead
# zeta(3,2) via known: sum_{n>m} 1/(n^3 m^2): converges polynomially! 900 terms ~ 1e-? tail ~ sum_{n>900} 1/n^3 * zeta(2) ~ 1e-6/2 — BAD. Use Euler: zeta(3,2) = 3*zeta(2)*zeta(3)/... no. Use the MPL-integral: zeta(3,2)=G(0,0,1,0,1;1)
# faster: sum to N=2000 then Richardson? Actually use the double integral directly is messy; use known value via mp.zeta? No mzeta in mpmath.
# Use the stuffle: zeta(3,2) + zeta(2,3) = z2*z3 - z5; zeta(2,3) = sum 1/(n^2 m^3) equally slow.
# Better: integral representation via my G-table is easiest:
print("Li32 time %.1fs" % (time.time()-t0))
K = (-8)*mp.mpf(pickle.load(open("atom_vals.pkl","rb"))['L32_b1']) - 32*Li32(-b, mp.mpf(1)) - 16*Li32(b, -mp.mpf(1)) - 16*Li32(-b, -mp.mpf(1)) + Li32(-b*b, mp.mpf(1))
print("K (symbol-kernel Li32 combo) =", mp.nstr(K, 20))
pickle.dump({'K': mp.nstr(K,240)}, open("K_val.pkl","wb"))
