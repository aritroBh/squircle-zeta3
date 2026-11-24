import sys, os, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
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

PHI = (Fr(1,2), Fr(1,2)); TWO = (Fr(2), Fr(0))
t0 = time.time()
# golden corner c0 = sqrt5 - 2 = phi^-3
P3 = W5Point("c=sqrt5-2 (golden3)", [-5, 0, 1], (Fr(-2), Fr(1)), basis=[TWO, PHI])
P3.S0 = [2]; P3.dimV = 2
lv3 = P3.letter_vectors()
print("letters:", {k: [str(x) for x in v] for k, v in lv3.items()}, flush=True)
S3 = P3.omega5()
SE3 = s01.ft_expand_sym(S3, 2)
obs3 = s01.obstruction_fixed(SE3)
print(f"obs3 [Q, T] = {[[str(x) for x in oo] for oo in obs3]}  ({time.time()-t0:.0f}s)", flush=True)

# controls: square and hexagon
P4 = W5Point("c=1/sqrt2 (square)", [-2, 0, 1], (Fr(0), Fr(1, 2)),
             basis=[(Fr(2), Fr(0)), (Fr(1), Fr(1))])
obs4 = s01.obstruction_fixed(s01.ft_expand_sym(P4.omega5(), 2))
print("square control:", [[str(x) for x in oo] for oo in obs4], "target [(84,0),(0,-38)]", flush=True)
P6 = W5Point("c=1/2 (hexagon)", [-1, 1], (Fr(1, 2),))
obs6 = s01.obstruction_fixed(s01.ft_expand_sym(P6.omega5(), 2))
print("hexagon control:", [[str(x) for x in oo] for oo in obs6], "target [(41,),(-33/2,)]", flush=True)
print(f"({time.time()-t0:.0f}s)")
