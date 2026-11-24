import time, sys
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L

mp.mp.dps = 1100
c3 = mp.mpf(1)/3
t0 = time.time()
M3q = L.corner_coeffs_hp(2*mp.acos(c3), 3)[3]
t1 = time.time()
print(f"quadrature dps 1100: {t1-t0:.0f}s")
open("c13_M3_quad_1100.txt", "w").write(mp.nstr(M3q, 1095))
M3ref = mp.mpf(open("c13_M3_hpl_1150.txt").read())
print("|quad(1100) - HPL(1150)| =", mp.nstr(abs(M3q - M3ref), 3))
