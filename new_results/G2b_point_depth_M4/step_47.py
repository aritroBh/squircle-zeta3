
import sys, os
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
_spec = importlib.util.spec_from_file_location("g2b_step05", G2B + "step_05.py")
s05 = importlib.util.module_from_spec(_spec)
s05.__name__ = "g2b_step05"
_spec.loader.exec_module(s05)
C0 = Fr(1,2)
w = (1,0,1,0,1)
Delta1 = s05.delta_word_sym(w, 0, C0)
S5w = s05.S_J(w, C0, 0)
D5w = s05.delta_minus_5sym(S5w)
print("delta_word_sym:", len(Delta1), " delta^-(5sym):", len(D5w))
keys = sorted(set(Delta1) | set(D5w))
for k in keys:
    a = Delta1.get(k, Fr(0)); b = D5w.get(k, Fr(0))
    ratio = (a / b) if b else None
    print(f"  {k}: delta_word={a}  delta^-={b}  ratio={ratio}")
