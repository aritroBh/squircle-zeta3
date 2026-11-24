import os, time, json
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
print('w12 planes progress:', len(json.load(open('w12_planes.json'))) if os.path.exists('w12_planes.json') else 0)
t0=time.time()
from fractions import Fraction as F
import itertools
# linear forms in (s,t): dict {(ca,cb,const)} with Fraction coeffs; represent as tuple (ca,cb,const)
def Ladd(p,q): return tuple(x+y for x,y in zip(p,q))
def Lsub(p,q): return tuple(x-y for x,y in zip(p,q))
def Leq(p,q): return all(x==y for x,y in zip(p,q))
ZERO=(F(0),F(0),F(0)); ONE=(F(0),F(0),F(1))
def T1(P):
    a,b,c,e,f=P
    sig=Lsub(Ladd(e,f),Ladd(Ladd(a,b),c))
    # 3F2(e-a, f-a, sig; sig+b, sig+c; 1)
    return (Lsub(e,a), Lsub(f,a), sig, Ladd(sig,b), Ladd(sig,c))
def perms(P):
    a,b,c,e,f=P
    out=[]
    for (x,y,z) in set(itertools.permutations((a,b,c))):
        for (u,v) in [(e,f),(f,e)]:
            out.append((x,y,z,u,v))
    return out
def orbit(P0):
    seen={P0}
    frontier=[P0]
    while frontier:
        P=frontier.pop()
        for Q in perms(P)+[T1(P)]:
            # normalize: tuples hashable
            if Q not in seen:
                seen.add(Q); frontier.append(Q)
    return seen
def conditions(P):
    a,b,c,e,f=P
    res=[]
    # collapse: numerator param zero
    if a==ZERO or b==ZERO or c==ZERO: res.append('TRIVIAL')
    # cancellation -> 2F1
    if Leq(a,e) or Leq(a,f) or Leq(b,e) or Leq(b,f) or Leq(c,e) or Leq(c,f): res.append('CANCELS->2F1')
    sig=Lsub(Ladd(e,f),Ladd(Ladd(a,b),c))
    if sig==ONE: res.append('SAALSCHUTZ')
    # Dixon: exists main A among {a,b,c}: denom = 1+A-B, 1+A-C
    for A,B,C in [(a,b,c),(b,a,c),(c,a,b)]:
        if (Leq(e,Lsub(Ladd(ONE,A),B)) and Leq(f,Lsub(Ladd(ONE,A),C))) or (Leq(f,Lsub(Ladd(ONE,A),B)) and Leq(e,Lsub(Ladd(ONE,A),C))):
            res.append('DIXON')
    # Watson: (1+A+B)/2 and 2C
    for A,B,C in [(a,b,c),(b,a,c),(c,a,b),(a,c,b),(b,c,a),(c,b,a)]:
        h=Ladd(ONE,Ladd(A,B)); h=(h[0]/2,h[1]/2,h[2]/2)
        d2=(C[0]*2,C[1]*2,C[2]*2)
        if (Leq(e,h) and Leq(f,d2)) or (Leq(f,h) and Leq(e,d2)):
            res.append('WATSON')
    # Whipple: {A,1-A} in numerator, e, 1+2C-e in denom for some C
    for A in (a,b,c):
        oneA=Lsub(ONE,A)
        if any(Leq(oneA,X) for X in (a,b,c)):
            for Cc in (a,b,c):
                for E in (e,f):
                    other=f if E is e else e
                    # need 1+2Cc-E == other
                    t2=Ladd(ONE,(Cc[0]*2,Cc[1]*2,Cc[2]*2))
                    if Leq(Lsub(t2,E),other):
                        res.append('WHIPPLE')
    return res
loci={'Pi1':((F(0),F(-1),F(0)), (F(0),F(1),F(0)), (F(1),F(1),F(1)), (F(0),F(1),F(1)), (F(2),F(1),F(1))),
      'Pi2':((F(0),F(-1),F(0)), (F(1),F(0),F(0)), (F(2),F(0),F(1)), (F(1),F(0),F(1)), (F(2),F(1),F(1))),
      'diag':((F(0),F(-1),F(0)), (F(0),F(1),F(0)), (F(2),F(0),F(1)), (F(1),F(0),F(1)), (F(3),F(0),F(1)))}
# NOTE: diag uses s=t=u: params (-u,u,1+2u;1+u,1+3u) : (a,b,c,e,f) = ((0,-1,0),(0,1,0),(2,0,1),(1,0,1),(3,0,1))
for name,P0 in loci.items():
    orb=orbit(P0)
    hits=set()
    for P in orb:
        for c in conditions(P):
            hits.add(c)
    print('%s: orbit size %d, summability hits: %s'%(name,len(orb),hits if hits else 'NONE'))
print('t=%.1fs'%(time.time()-t0))
