import sys, pickle
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
d=pickle.load(open(WORK+'/cands29_reduced.pickle','rb'))
fw=d['freewords']; red=d['cands']
from flint import nmod_mat
import random

def solve_membership(p, rpk, rpkl):
    R=pickle.load(open(BASE+f'/research/08_theory/A3h_weight4_target/rref_T3_{p}.pkl','rb'))
    r=R['r'] if rpk else (-R['r'])%p
    nfD=[(R['nf']['Drat'][j]+r*R['nf']['Dsq'][j])%p for j in range(110)]
    nfP=R['nf']['P4']
    # matrix A: 110 x 30 : [nfP, N_1..N_29]
    # choose 30 independent rows
    Acols=[[nfP[j]%p for j in range(110)]]+[[ (red[k][j]%p) for j in range(110)] for k in range(29)]
    # find 30 rows where invertible: greedy pivot
    M=nmod_mat(110,30,p)
    for j in range(110):
        for c in range(30): M[j,c]=Acols[c][j]
    # solve M x = nfD: use rref on augmented
    Mg=nmod_mat(110,31,p)
    for j in range(110):
        for c in range(30): Mg[j,c]=Acols[c][j]
        Mg[j,30]=nfD[j]
    RR,rank=Mg.rref(inplace=True)
    # check consistency: rows >=30 must be zero in first 30 cols -> if rank>30 inconsistent
    if rank>30:
        return None, rank
    # solution: pivot rows
    x=[0]*30
    for i in range(rank):
        # find pivot col
        c=0
        while not RR[i,c]: c+=1
        if c<30: x[c]=int(RR[i,30])
    return x, rank

for p in (1000000009,1000000033):
    for rpk in (True,False):
        x,rank=solve_membership(p,rpk,None)
        R=pickle.load(open(BASE+f'/research/08_theory/A3h_weight4_target/rref_T3_{p}.pkl','rb'))
        rused=R['r'] if rpk else (-R['r'])%p
        if x is None:
            print(f'p={p} r={rused}: INCONSISTENT (rank {rank})')
        else:
            c=x[0]; y=x[1:]
            print(f'p={p} r={rused}: c={c}, y[:6]={y[:6]}, max|y| residues={max(yy if yy<p//2 else p-yy for yy in y)}')
