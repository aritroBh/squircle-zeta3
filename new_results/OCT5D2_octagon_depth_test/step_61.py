import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, add, scale, obstruction

data = pickle.load(open("symbol_fit2.pkl","rb"))
words, M, b, cols = data['words'], data['M'], data['b'], data['cols']
x = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
R2 = Q(0,1); B = Q(-1,1)
def w32(xx, yy): return (ZERO, ZERO, ONE/xx, ZERO, ONE/(xx*yy))
name2qs = {'b':B,'-b':-B,'1/sqrt2':R2/2,'-b^2':-B*B,'2-sqrt2':2-R2,'1-1/sqrt2':ONE-R2/2}
beta = [(label, xv) for (label, kind), xv in zip(cols, x) if kind=='li32' and xv != ZERO]
print("beta part:", [(l, str(v)) for l,v in beta])
sbeta = {}
for label, xv in beta:
    xs, ys = label[5:].split(',')
    xx = name2qs[xs]; yy = ONE if ys=='1' else -ONE
    sbeta = add(sbeta, scale(S.gsym(w32(xx,yy)), xv))
obs_beta = obstruction(sbeta)
print("pair of beta (Li32) part:", obs_beta)
print("expected obstruction:    (-1485+1035*sqrt2, 900-540*sqrt2)")
print("match:", obs_beta == [Q(-1485,1035), Q(900,-540)])

# also verify the FULL solution symbol equals W5 total symbol (rebuild column symbols quickly)
# do it via matrix: M x == b exactly (already by construction) - just re-verify explicitly
Mxv = [sum((M[i][c]*x[c] for c in range(len(x))), ZERO) for i in range(len(words))]
print("M x == b exactly:", all(Mxv[i] == b[i] for i in range(len(words))))
# classical part pair must be zero: compute pair of (W5sym - sbeta) = classical symbol
raw = pickle.load(open("block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
W5sym = parse(raw['total'])
classical_sym = add(W5sym, scale(sbeta, -1))
print("pair of classical part:", obstruction(classical_sym))
