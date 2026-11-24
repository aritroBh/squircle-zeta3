import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, add, scale, obstruction

data = pickle.load(open("symbol_fit2.pkl","rb"))
words, M, b, cols = data['words'], data['M'], data['b'], data['cols']
x = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']

# reconstruct: full solution symbol, classical part, li32 part
# need column symbols again: reload from symbol_fit2 construction -> we saved only labels; reconstruct via pickle of columns? we didn't save symbols.
# quick: recompute li32 pair of the beta part from the labels
from fractions import Fraction as F
R2 = Q(0,1); B = Q(-1,1)
def w32(xx, yy): return (ZERO, ZERO, ONE/xx, ZERO, ONE/(xx*yy))
LI32 = {'(b,1)':(B,ONE),'(1/sqrt2,1)':(R2/2,ONE),'(-b,1)':(-B,ONE),
        '(b,-1)':(B,-ONE),'(-b,-1)':(-B,-ONE),'(1/sqrt2,-1)':(R2/2,-ONE),
        '(1-1/sqrt2,1)':(ONE-R2/2,ONE),'(2-sqrt2,1)':(2-R2,ONE),'(-b^2,1)',:None}
# simpler: collect beta from solution
beta = [(label, xv) for (label, kind), xv in zip(cols, x) if kind=='li32' and xv != ZERO]
print("beta part:", [(l, str(v)) for l,v in beta])
sbeta = {}
for label, xv in beta:
    import ast
    # parse (x,y) from label Li32(b,1) style: map names
    arg = label[5:]
    name2qs = {'b':B,'-b':-B,'1/sqrt2':R2/2,'-b^2':-B*B,'2-sqrt2':2-R2,'1-1/sqrt2':ONE-R2/2}
    xs, ys = arg.split(',')
    xx = name2qs[xs]; yy = ONE if ys=='1' else -ONE
    sbeta = add(sbeta, scale(S.gsym(w32(xx,yy)), xv))
obs_beta = obstruction(sbeta)
print("pair of beta (Li32) part:", obs_beta)
print("expected obstruction: (-1485+1035*sqrt2, 900-540*sqrt2)")
print("match:", obs_beta == [Q(-1485,1035), Q(900,-540)])
