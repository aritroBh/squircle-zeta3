import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
open('progress.log','a').write('step_76-81: diag vanishing W10/11 confirmed; values at (3,1,2,5),(1,1,1,2),(1,2,2,1); diag c12 coords (z39=41328,z57=-26568 ratio -14/9, w*=0); diag c12 Cauchy |diff|=4.84e-71 PASS\n')
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, itertools
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w=sp.symbols('x y z w')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
Q4={tuple(xx-1 for xx in e):c for e,c in forms[(8,3)].items()}
C3=sp.cancel(10*to_sp(Q4)/(x+y+z+w))
# invariant cubics: orbit sums under sigma=(x w)(y z)
def orb(e):
    e2=e[::-1]
    return x**e[0]*y**e[1]*z**e[2]*w**e[3] + (0 if e2==e else x**e2[0]*y**e2[1]*z**e2[2]*w**e2[3])
orbs={}
for e in itertools.product(range(4), repeat=4):
    if sum(e)==3:
        k=tuple(sorted([e,e[::-1]]))
        if k not in orbs: orbs[k]=orb(e)
basis=list(orbs.items())
print('invariant cubics:', len(basis))
# solve C3 in this basis
M=[]
for k,p in basis:
    d={}
    for mon,c in sp.Poly(p).terms(): d[mon]=c
    M.append(d)
allmons=sorted(set(m for d in M for m in d))
Amat=[[d.get(m,0) for d in M] for m in allmons]
C3d={mon:c for mon,c in sp.Poly(C3).terms()}
b=[C3d.get(m,0) for m in allmons]
sol=sp.linsolve((sp.Matrix(Amat), sp.Matrix(b)))
sol=list(sol)[0]
print('C3 in invariant basis:')
for (k,p),c in zip(basis,sol):
    if c!=0: print('  %s * %s  where orbit=%s'%(c,k,p))
print('t=%.1fs'%(time.time()-t0))
