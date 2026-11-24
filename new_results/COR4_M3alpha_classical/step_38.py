import sys, json, re
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R2 = json.load(open(G1 + "hpl_fit_k2_(0,0,0).json"))
print(len(R2), R2[:8])
