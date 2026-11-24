
import sys
sys.path.insert(0, '<project>/research/08_theory/COR8_depth_two_any_angle')
import mpmath as mp
mp.mp.dps = 50
import sympy as sp
from cor8_series import *
c = sp.Symbol('c')
S1 = sscale(slog_of_rat(1 - (1/(1-c)), c), -1)
print("S1:", sorted(S1.items()))
dlogg = sdiff(slog_of_rat(1/(1-c), c))
print("dlogg:", sorted(dlogg.items())[:6])
P = smul(S1, dlogg)
print("P:", sorted(P.items())[:8])
I = sintegrate(P)
print("I:", sorted(I.items())[:8])
