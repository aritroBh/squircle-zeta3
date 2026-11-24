import sys, os, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
GOLD = BASE + "/GOLD5_golden_corner_w5/"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); mod.__name__ = name
    spec.loader.exec_module(mod)
    return mod
s01 = load("g2b_step01", G2B + "step_01.py")
from step_13 import W5Point

t0 = time.time()
# golden corner c0 = 1/phi = (-1+sqrt5)/2, basis {2, phi}
PHI = (Fr(1,2), Fr(1,2))
TWO = (Fr(2), Fr(0))
P = W5Point("c=1/phi (golden1)", [-5, 0, 1], (Fr(-1,2), Fr(1,2)), basis=[TWO, PHI])
P.S0 = [2]; P.dimV = 2
print("field deg", P.F.d, "dimV patched", P.dimV)
# check letter vectors
lv = P.letter_vectors()
print("letters:", {k: [str(x) for x in v] for k, v in lv.items()})
# expect: c=phi^-1 -> (0,-1); 1-c=phi^-2 -> (0,-2); 1+c=phi -> (0,1)
print(f"({time.time()-t0:.0f}s)")
