# chkA_tower5.py -- own check of hypothesis of Lemma 3 at n=5: tower law for keys (5,0),(4,1),(3,2),
# depths 0..6, from two independent exact data sources (TOWN1 truncated engine k<=52; WTS full reducer towers k<=14).
import pickle, sys, os
from fractions import Fraction as Fr
R='<project>/research/08_theory/'
T1=pickle.load(open(R+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
WT={}
for k in range(2,15):
    WT[k]=pickle.load(open(R+f'WTS_weight_k_structure/theory/tower_{k}.pkl','rb'))
def fd(seq,m):
    s=list(seq)
    for _ in range(m): s=[s[i+1]-s[i] for i in range(len(s)-1)]
    return s
def ok_poly(seqs, deg):
    """seqs: list of 4 component sequences (consecutive k). returns (is_poly_deg<=deg, nchecks, leading_nonzero)"""
    nchk=len(seqs[0])-(deg+1)
    allzero=all(all(x==0 for x in fd(s,deg+1)) for s in seqs)
    lead=any(fd(s,deg)[0]!=0 for s in seqs) if len(seqs[0])>deg else None
    # exact degree: deg-th difference nonzero (constant) -> check all equal & nonzero
    return allzero, nchk, lead
n=5
keys=[(5,0),(4,1),(3,2)]
bad=0; tot=0
for key in keys:
    for delta in range(0,7):
        P=delta//2
        for k0 in (5,6):
            # TOWN1
            ks=list(range(k0,53))
            seqs=[[T1[key][delta][k][c] for k in ks] for c in range(4)]
            ap,nchk,lead=ok_poly(seqs,P)
            tot+=1
            if not (ap and lead):
                bad+=1; print('BAD TOWN1',key,delta,k0,ap,lead)
            # WTS full-reducer towers
            r=2*n-2-delta
            ks2=list(range(k0,15))
            seqs2=[[WT[k][key].get(r,(Fr(0),)*4)[c] for k in ks2] for c in range(4)]
            ap2,nchk2,lead2=ok_poly(seqs2,P)
            tot+=1
            if not (ap2 and lead2):
                bad+=1; print('BAD WTS',key,delta,k0,ap2,lead2)
            # cross-source equality
            for k in ks2:
                a=tuple(T1[key][delta][k]); b=tuple(WT[k][key].get(r,(Fr(0),)*4))
                if a!=b: bad+=1; print('MISMATCH',key,delta,k,a,b)
            if k0==6: print(f'{key} delta={delta} deg={P}: TOWN1 k=6..52 ok={ap} (checks {nchk}), WTS k=6..14 ok={ap2} (checks {nchk2}), exact-degree={lead and lead2}')
print('total tests',tot,'bad',bad)
# Does k=5 lie on the polynomial fitted from k>=6? (T claims k>=n=5 for off-diagonal): test explicitly via extrapolation
def extrap(vals,ks,k):
    # Lagrange through points (ks[i],vals[i]) evaluated at k, exact
    tot=Fr(0)
    for i,xi in enumerate(ks):
        li=Fr(1)
        for j,xj in enumerate(ks):
            if i!=j: li*=Fr(k-xj,xi-xj)
        tot+=vals[i]*li
    return tot
offs=0
for key in keys:
    for delta in range(0,7):
        P=delta//2
        ks=list(range(6,6+P+1))
        for c in range(4):
            vals=[T1[key][delta][k][c] for k in ks]
            if extrap(vals,ks,5)!=T1[key][delta][5][c]: offs+=1; print('k=5 off poly',key,delta,c)
print('k=5 off-polynomial components (fit from 6..):',offs)
# the diagonal parent (3,3): off at k=6, on from k>=7 (shift)
print('diag (3,3) check:')
for delta in range(0,9):
    P=delta//2
    if delta not in T1[(3,3)]: continue
    seq=lambda k0:[[T1[(3,3)][delta][k][c] for k in range(k0,53)] for c in range(4)]
    a6=ok_poly(seq(6),P)[0]; a7=ok_poly(seq(7),P)[0]
    print(f' (3,3) delta={delta}: polynomial from k=6? {a6}; from k=7? {a7}')
