
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 65
import reflib as R
# independent closed-form checks of the trailing-zero shuffle recursion
v1 = R.hpl_value((1,0), mp.mpf(1)/2)
f1 = -mp.log(2)**2 - R.hpl_value((0,1), mp.mpf(1)/2)   # H_{1,0} = log(1/2)(-log(1/2))... = -log^2 2 - Li_2(1/2)? log x*(-log(1-x)) - Li_2(x)
print("H_10:", mp.nstr(v1, 10), "expect", mp.nstr(f1, 10), "err", mp.nstr(abs(v1-f1), 3))
v2 = R.hpl_value((0,1,0), mp.mpf(1)/2)
f2 = (R.hpl_value((0,1), mp.mpf(1)/2)*mp.log(mp.mpf(1)/2) - R.hpl_value((0,0,1), mp.mpf(1)/2))/2
print("H_010:", mp.nstr(v2, 10), "expect", mp.nstr(f2, 10), "err", mp.nstr(abs(v2-f2), 3))
