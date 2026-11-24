# step_60: inspect W_parts_e0 structure
import os, sys, pickle
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9b_octagon_debug"))
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9c_octagon_close"))
import oct9c_common as oc
import sympy as sp
D = oc.load_all()
W3e, Wloe = D['Wp0']
print(type(W3e), type(Wloe))
print("W3e:", repr(W3e)[:500])
print("Wloe:", repr(Wloe)[:500])
