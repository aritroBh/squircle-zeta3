import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_common.py").read()
src = src.replace(
    "    h1,h2,h3 = (Hder(m)/H0/mp.factorial(m) for m in range(1,4))\n",
    "    h1,h2,h3 = (omega**m * Hder(m)/H0/mp.factorial(m) for m in range(1,4))\n")
open("dcy2_common.py","w").write(src)

from mpmath import mp, mpf, pi, sqrt, exp
mp.dps = 60
import importlib, dcy2_common
importlib.reload(dcy2_common)
from dcy2_common import model, target, leading, I_tilde

for N in [100, 400, 2500, 25600, 102400]:
    N = mpf(N)
    Mo = model(N); T1 = target(N); L = leading(N)
    print(f"N={int(N):7d}  (model-T1)/(L/N) = {float(abs(Mo-T1)/(L/N)):8.4f}")
for N in [100, 400, 2500]:
    N = mpf(N)
    It = I_tilde(N); Mo = model(N); T1 = target(N); L = leading(N)
    print(f"N={int(N):7d}  |It-model|/(L/N) = {float(abs(It-Mo)/(L/N)):8.4f}  |It-T1|/(L/N) = {float(abs(It-T1)/(L/N)):8.4f}")
