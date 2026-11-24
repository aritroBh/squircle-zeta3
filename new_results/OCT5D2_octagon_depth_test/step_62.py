import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, add, scale, obstruction
from fractions import Fraction as F

R2 = Q(0,1); B = Q(-1,1)
def w32(xx, yy): return (ZERO, ZERO, ONE/xx, ZERO, ONE/(xx*yy))
for nm,(xx,yy) in {'(b,1)':(B,ONE),'(1/sqrt2,1)':(R2/2,ONE),'(-b,1)':(-B,ONE)}.items():
    print("pair Li32%s =" % nm, obstruction(S.gsym(w32(xx,yy))))

c1 = Q(-3600,2160); c2 = Q(11880,-8280)
sb = add(scale(S.gsym(w32(B,ONE)), c1), scale(S.gsym(w32(R2/2,ONE)), c2))
print("beta pair:", obstruction(sb))
# expected: c1*(0,-1/4) + c2*(-1/8,0)
print("manual: A =", c2*Q(F(-1,8)), " B =", c1*Q(F(-1,4)))
raw = pickle.load(open("block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
W5sym = parse(raw['total'])
print("W5 total pair:", obstruction(W5sym))
classical_sym = add(W5sym, scale(sb, -1))
print("classical-part pair (W5 - beta):", obstruction(classical_sym))
