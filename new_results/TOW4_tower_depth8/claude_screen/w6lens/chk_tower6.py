import sys, os, pickle
sys.path.insert(0, os.getcwd())
from k_load import *
from fractions import Fraction as Fr
from math import comb
H = pickle.load(open(TH+'/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def Wtower(k):
    return pickle.load(open(f'{TH}/WTS_weight_k_structure/theory/tower_{k}.pkl','rb'))
def lagr_vec(xs, ys, x):
    s=[Fr(0)]*4
    for i,xi in enumerate(xs):
        L=Fr(1)
        for j,xj in enumerate(xs):
            if i!=j: L*=Fr(x-xj, xi-xj)
        for t in range(4): s[t]+=ys[i][t]*L
    return tuple(s)
def fd(seq, n):
    out=[Fr(0)]*4
    for j in range(n+1):
        c=(-1)**(n-j)*comb(n,j)
        for t in range(4): out[t]+=c*seq[j][t]
    return tuple(out)
keys6=[(6,0),(5,1),(4,2),(3,3)]
print("== (1a) TOWN1 data: degree-6 keys, depth 0..8, k from k0 to 52")
allok=True
for key in keys6:
    n=sum(key); diag=(key[0]==key[1])
    k0 = n+1 if diag else n
    row=[]
    for dl in range(0,9):
        P=dl//2
        seq=[H[key][dl][k] for k in range(k0,53)]
        # (P+1)-th difference at every window vanish, P-th nonzero at some window
        vanish=True; npts=len(seq)
        for s in range(0,len(seq)-(P+1)):
            if any(c!=0 for c in fd(seq[s:s+P+2],P+1)): vanish=False
        nonzero=any(any(c!=0 for c in fd(seq[s:s+P+1],P)) for s in range(0,len(seq)-P))
        # leading coefficient: P-th difference constant = P! * lead
        row.append((dl,P,vanish,nonzero,npts-(P+1)))
        if not (vanish and nonzero): allok=False
    print(key,'k0=',k0, [(r[0],r[1],'OK' if (r[2] and r[3]) else 'FAIL', r[4]) for r in row])
print('all degree-6 (key,depth) cells: polynomial of degree exactly floor(delta/2) from k0..52 :', allok)

print("\n== (1b) diagonal (3,3): is k=6 on the stable polynomial (fit k=7..)?")
for dl in range(0,9):
    P=dl//2
    xs=list(range(7,7+P+1)); ys=[H[(3,3)][dl][k] for k in xs]
    pred=lagr_vec(xs,ys,6)
    act=H[(3,3)][dl][6]
    print('delta',dl,'deg',P,'k=6 on polynomial?', pred==act, ' disc=', tuple(str(a-b) for a,b in zip(act,pred)))
print("\n== (1c) off-diagonal degree-6: is k=6 ON the polynomial (fit from k=7..)?")
for key in [(6,0),(5,1),(4,2)]:
    res=[]
    for dl in range(0,9):
        P=dl//2
        xs=list(range(7,7+P+1)); ys=[H[key][dl][k] for k in xs]
        res.append(lagr_vec(xs,ys,6)==H[key][dl][6])
    print(key,res)

print("\n== (1d) independent full-reducer towers (WTS tower_k.pkl, k<=14) vs TOWN1 for degree-6 keys, delta 0..8")
tot=0;bad=0
for k in range(6,15):
    T=Wtower(k)
    for key in keys6:
        if key not in T: print('missing',k,key); continue
        for dl in range(0,9):
            r=2*6-2-dl
            got=T[key].get(r,(Fr(0),)*4)
            want=H[key][dl][k]
            tot+=1
            if tuple(got)!=tuple(want): bad+=1
print('cells',tot,'bad',bad)
# degree-of-k check using ONLY the independent WTS data (k<=14): fit from first P+1, check rest
print("\n== (1e) same law from WTS-only data (k=6/7..14)")
ok=True
for key in keys6:
    diag=(key[0]==key[1]); k0=7 if diag else 6
    for dl in range(0,9):
        P=dl//2
        seq=[Wtower(k)[key].get(10-dl,(Fr(0),)*4) for k in range(k0,15)]
        if len(seq)<P+3: print('too few points',key,dl); continue
        v=all(not any(c!=0 for c in fd(seq[s:s+P+2],P+1)) for s in range(len(seq)-P-1))
        nz=any(any(c!=0 for c in fd(seq[s:s+P+1],P)) for s in range(len(seq)-P))
        if not(v and nz): ok=False; print('FAIL',key,dl)
print('WTS-only degree-6 tower law OK:',ok)
