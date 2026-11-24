import os, sys, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, const_pow, obstruction, li_symbol, gsym
from fractions import Fraction as F

R2 = Q(0,1); B = Q(-1,1)
# S-unit arguments where z and 1-z are both S-units (test via li_symbol working)
cands = {'b': B, '-b': -B, 'b^2': B*B, '-b^2': -B*B, '1/2': Q(F(1,2)), '1/sqrt2': R2/2,
         'sqrt2': R2, '-sqrt2': -R2, '2-sqrt2': 2-R2, 'p=2+sqrt2': P_LET, 'q=-(1+sqrt2)': Q_LET,
         '-1': -ONE, '2': Q(2), '1/4': Q(F(1,4)), '4': Q(4), '1/(2 sqrt2)': R2/4,
         '2 sqrt2': 2*R2, '1/(1+sqrt2)': ONE/(1+R2), '-(2-sqrt2)': -(2-R2), 'sqrt2-1 over 2?': (R2-1)/2}
zs = []
for name, z in cands.items():
    try:
        li_symbol(2, z); zs.append((name, z))
    except Exception as e:
        print("skip", name, ":", e)
print("usable Li_n arguments:", [n for n,_ in zs])

# quick check of symbol of W5 total from pickle
raw = pickle.load(open("block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
W5sym = parse(raw['total'])
print("W5 total symbol words:", len(W5sym))
