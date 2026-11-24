import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, symengine as S
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
P8=forms[(8,3)]
Q4={tuple(x-1 for x in e):c for e,c in P8.items()}
# reversal invariance: (e1,e2,e3,e4)->(e4,e3,e2,e1)
rev = all(Q4.get(e[::-1])==c for e,c in Q4.items())
print('Q4 reversal invariant:', rev)
# convert to sympy
import sympy as sp
x,y,z,w=sp.symbols('x y z w')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
q4=to_sp(Q4)
print('factoring...')
t=time.time()
fac=sp.factor(q4)
print('factor(q4) =', fac, ' (%.0fs)'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
