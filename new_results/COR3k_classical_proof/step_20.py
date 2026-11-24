import sympy as sp
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
print("loaded", len(CD), "entries")
print("C[0,0,1] =", CD['[0, 0, 1]'])
print("C[0,1,-1] has", len(CD['[0, 1, -1]'].args), "additive terms")

# test dpoly on a Li3 atom
e = sp.polylog(3, ARGS['(1-x)/(1+x)'])
d = dpoly(e)
print("d/dx Li3((1-x)/(1+x)) =", sp.cancel(d))
# expect: polylog(2, arg) * arg'/arg
chk = sp.cancel(d - sp.polylog(2, ARGS['(1-x)/(1+x)'])*sp.diff(ARGS['(1-x)/(1+x)'], x)/ARGS['(1-x)/(1+x)'])
print("matches chain rule:", chk == 0)
