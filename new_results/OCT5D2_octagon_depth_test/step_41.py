import os, sys, pickle, json
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, obstruction, gsym
from fractions import Fraction as F

R2 = Q(0,1); B = Q(-1,1)   # sqrt2, b = sqrt2-1
def w32(x, y): return (ZERO, ZERO, ONE/x, ZERO, ONE/(x*y))
atoms = {'(b,1)': (B, ONE), '(1/sqrt2,1)': (R2/2, ONE), '(-b,1)': (-B, ONE), '(b,-1)': (B, -ONE),
         '(1/sqrt2,-1)': (R2/2, -ONE), '(1-1/sqrt2,1)': (ONE-R2/2, ONE), '(2-sqrt2,1)': (2-R2, ONE),
         '(-b^2,1)': (-B*B, ONE)}
pairs = {}
for name,(x,y) in atoms.items():
    s = gsym(w32(x,y))
    obs = obstruction(s)
    pairs[name] = obs
    print(f"Li32{name}: symbol words {len(s):3d}  pair = ({obs[0]}, {obs[1]})")

# N5 combination check: must pair to (84, -38 sqrt2)
N5 = [('(b,1)', Q(F(779), F(-695,56))), ('(1/sqrt2,1)', Q(F(-672))), ('(-b,1)', Q(F(-374), F(76))), ('(b,-1)', Q(F(-405), F(4951,56)))]
tot = {}
for name, c in N5:
    tot = add(tot, scale(gsym(w32(*atoms[name])), c))
print("N5 pair:", obstruction(tot), " (expect 84, -38*sqrt2)")
