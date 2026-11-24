
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "g2b_dbg.py", "w").write('''
import sys, os
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    mod.__name__ = name
    spec.loader.exec_module(mod)
    return mod
s05 = load("g2b_step05", G2B + "step_05.py")
''')
print("helper written")
