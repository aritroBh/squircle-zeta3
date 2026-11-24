
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
import types
# prevent __main__ block: temporarily set __name__
s05.__name__ = "g2b_step05"
_spec.loader.exec_module(s05)

C0 = Fr(1,2)
# single word: w = (0,0,0,0,1) with coefficient 1
w = (0,0,0,0,1)
Delta1 = s05.delta_word_sym(w, 0, C0)
print("delta_word_sym of J((0,0,0,0,1); 0->1/2):")
for k, v in sorted(Delta1.items()):
    print("   ", k, v)
# 5-symbol of the same word and its delta^-
S5w = s05.S_J(w, C0, 0)
D5w = s05.delta_minus_5sym(S5w)
print("delta^- of S(J((0,0,0,0,1);1/2)):")
for k, v in sorted(D5w.items()):
    print("   ", k, v)
