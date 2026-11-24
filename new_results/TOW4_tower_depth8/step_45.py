import sys
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
d, m = T4.d, T4.m
for delta in range(5,9):
    for s in range(0, delta//2+1):
        f = T4.forced(delta,s)
        num, den = sp.fraction(sp.cancel(f))
        nd = sp.Poly(num, d).degree() if num!=0 else -1
        nt = len(sp.Add.make_args(sp.expand(num)))
        print(f"F({delta},{s}): num deg_d={nd} terms={nt} den={sp.factor(den)}")
