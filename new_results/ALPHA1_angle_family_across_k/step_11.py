"""step_11.py -- weight two: T(alpha) = int_r^1 R_w2 lt dw, numeric + closed-form guess;
weight three: slope kernels S_20, S_11 (dK/dk) residues, and the slope value."""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50

def qpart(t):
    return Fr(*t[0]), Fr(*t[1])

# T(alpha) from the universal kernel
def T_alpha(cnum, cden):
    c = mp.mpf(cnum)/cden
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

# squircle control: c=1/sqrt2 -> T = (pi^2 + 6 log^2 2)/24
Tsq = T_alpha(1, mp.sqrt(2))
lam = mp.log(2)
print("squircle T =", mp.nstr(Tsq, 30), " target:", mp.nstr((mp.pi**2+6*lam**2)/24, 30),
      " diff:", mp.nstr(abs(Tsq-(mp.pi**2+6*lam**2)/24), 3))

for cnum, cden in ((1,2),(3,5),(1,3),(4,5)):
    T = T_alpha(cnum, cden)
    print(f"c={cnum}/{cden}: T = {mp.nstr(T, 40)}")

# guess: T = (pi^2 + A * log^2(something))/24?  Try T = (pi^2 + 6 L^2)/24 with L = log(2(1-r))? = log(4c/(1+c))
print()
for cnum, cden in ((1,2),(3,5),(1,3),(4,5)):
    c = mp.mpf(cnum)/cden; r = (1-c)/(1+c)
    T = T_alpha(cnum, cden)
    L = mp.log(2*(1-r))   # = log(4c/(1+c))
    cand = (mp.pi**2 + 6*L**2)/24
    print(f"c={cnum}/{cden}: T={mp.nstr(T,25)}  (pi^2+6 L^2)/24={mp.nstr(cand,25)} diff={mp.nstr(abs(T-cand),3)}")
print("STEP11 DONE")
