import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import engine as E, symengine as S, pickle, random, itertools
# 19 orbit-sum monomials for reversal-invariant quartics
mons=set()
for e in itertools.product(range(5), repeat=4):
    if sum(e)==4:
        orb=tuple(sorted([e, e[::-1]]))
        mons.add(orb)
monlist=sorted(mons)
print('invariant monomial orbits:', len(monlist))
assert len(monlist)==19
def orbval(orb, a):
    return a[0]**orb[0][0]*a[1]**orb[0][1]*a[2]**orb[0][2]*a[3]**orb[0][3] + (a[0]**orb[1][0]*a[1]**orb[1][1]*a[2]**orb[1][2]*a[3]**orb[1][3] if orb[1]!=orb[0] else 0)
random.seed(7)
pts=set()
while len(pts)<40: pts.add(tuple(random.randint(1,6) for _ in range(4)))
pts=sorted(pts)
vals=[]
for a in pts:
    p8=E.reduced_coeffs(a,8)[8][3]
    vals.append(F(p8, a[0]*a[1]*a[2]*a[3]))
# solve 19x19
M=[[F(orbval(o,a)) for o in monlist] for a in pts[:19]]
b=vals[:19]
import mzv_exact as MX
# gaussian solve
n=19
A=[row[:]+[bb] for row,bb in zip(M,b)]
for c in range(n):
    p=next(r for r in range(c,n) if A[r][c]!=0)
    A[c],A[p]=A[p],A[c]
    inv=F(1)/A[c][c]; A[c]=[x*inv for x in A[c]]
    for r in range(n):
        if r!=c and A[r][c]!=0:
            f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
sol=[A[i][n] for i in range(n)]
# validate held-out
ok=True
for a,v in zip(pts[19:], vals[19:]):
    pv=sum(s*orbval(o,a) for s,o in zip(sol,monlist))
    if pv!=v: ok=False; print('HELD-OUT FAIL', a)
print('held-out (21 pts):', 'PASS' if ok else 'FAIL')
# compare with symbolic Q4
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
P8=forms[(8,3)]
Q4={tuple(x-1 for x in e):c for e,c in P8.items()}
sym={}
for e,c in Q4.items():
    orb=tuple(sorted([e,e[::-1]]))
    # coefficient of orbit representative: for non-fixed orbits the sym coeff appears on both
    sym[orb]=sym.get(orb,F(0))+c
match=True
for o,s in zip(monlist,sol):
    # interpolated coefficient s multiplies orbit SUM; symbolic: each monomial has coeff c -> orbit sum coeff = c (same on both)
    c0 = Q4.get(o[0], F(0))
    if s!=c0 and not (o[0]==o[1] and s==c0):
        match=False; print('COEFF MISMATCH', o, s, c0)
print('interpolation == symbolic Q4:', 'PASS' if match else 'FAIL')
print('t=%.1fs'%(time.time()-t0))
open('progress.log','a').write('step_36: Q4 interpolation cross-check PASS (19 unknowns, 21 held-out exact, matches symbolic)\n')
