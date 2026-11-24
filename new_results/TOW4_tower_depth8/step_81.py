import sys, time
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
d, m = T4.d, T4.m
for dp in range(0,9):
    U = sp.cancel(sum((-1)**s * T4.forced(dp, s) for s in range(0, dp//2+1)))
    num, den = sp.fraction(U)
    print(f"U({dp}) = ( {sp.expand(num)} ) / ( {sp.factor(den)} )")
