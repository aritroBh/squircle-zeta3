import os, sys
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp
mp.dps = 60
from dcy2_common import model, target, leading

for N in [400, 1600, 6400, 25600, 102400]:
    N = mpf(N)
    Mo = model(N); T1 = target(N); L = leading(N)
    print(f"N={int(N):7d}  (model-T1)/(L/N) = {float(abs(Mo-T1)/(L/N)):10.4f}   "
          f"(model-T1)/(L/sqrtN) = {float(abs(Mo-T1)/(L/sqrt(N))):.5f}")
