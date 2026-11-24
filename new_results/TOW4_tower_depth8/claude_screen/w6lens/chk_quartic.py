import sys, os, json
sys.path.insert(0, os.getcwd())
from k_load import *
from math import comb
from fractions import Fraction as Fr

def vec(k, key, pole): return comp(k, key, pole)
def fdiff(seq, order):
    """forward difference (componentwise tuples) at start of seq"""
    n = order
    out = [Fr(0)]*4
    for j in range(n+1):
        c = (-1)**(n-j)*comb(n, j)
        for t in range(4): out[t] += c*seq[j][t]
    return tuple(out)

ks = list(range(5, 22))
# union of supports over k=5..21 and keys
allpoles = {}
for key in KEYS5:
    s=set()
    for k in ks: s |= poleset(k, key)
    allpoles[key]=sorted(s)
    print(key, 'pole union k=5..21:', [(str(a),str(b)) for a,b in allpoles[key]])
    for k in ks:
        ps = poleset(k,key)
        if ps != set(allpoles[key]) and k>=6: print('   support differs at k=',k, [(str(a),str(b)) for a,b in sorted(ps)])

# also: any component with imaginary part nonzero?
imag_nonzero = 0
for key in KEYS5:
    for k in ks:
        for p in allpoles[key]:
            v = vec(k,key,p)
            if v[2]!=0 or v[3]!=0: imag_nonzero += 1
print('nonzero imaginary components (real residues should be real):', imag_nonzero)

print("\n== (3a) 5-term identity K6 = 5K7 - 10K8 + 10K9 - 5K10 + K11, all poles of union support, all 4 components")
tot=0; bad=0
for key in KEYS5:
    for p in allpoles[key]:
        seq = [vec(k,key,p) for k in range(6,12)]
        d5 = fdiff(seq,5)  # K6 -5K7+10K8-10K9+5K10-K11  (= Delta^5 K6, sign (+)(-1)^{5-j})
        for t in range(4):
            tot+=1
            # Delta^5 = sum (-1)^{5-j} C(5,j) K_{6+j} = -K6+5K7-10K8+10K9-5K10+K11 ; zero iff identity
            if d5[t]!=0: bad+=1
        # explicit identity as stated
        for t in range(4):
            lhs = seq[0][t]; rhs = 5*seq[1][t]-10*seq[2][t]+10*seq[3][t]-5*seq[4][t]+seq[5][t]
            assert (lhs==rhs)==(d5[t]==0)
print('components', tot, 'bad', bad)

print("\n== (3b) Delta^5 on all windows k..k+5 for k=7..16 (so quartic through 7..11 reproduces 12..21)")
tot=0; bad=0
for key in KEYS5:
    for p in allpoles[key]:
        for k0 in range(7,17):
            seq=[vec(k,key,p) for k in range(k0,k0+6)]
            d5=fdiff(seq,5)
            for t in range(4):
                tot+=1; bad+= (d5[t]!=0)
print('window-components', tot, 'bad', bad)

# direct Lagrange prediction from 7..11 to 12..21 and to 6, 5
def lagr(xs, ys, x):
    s=Fr(0)
    for i,xi in enumerate(xs):
        L=Fr(1)
        for j,xj in enumerate(xs):
            if i!=j: L*= Fr(x-xj, xi-xj)
        s+= ys[i]*L
    return s
print("\n== (3c) Lagrange quartic through k=7..11 evaluated at 5,6,12..21, vs recorded")
res={}
for kk in [5,6]+list(range(12,22)):
    tot=0; bad=0
    for key in KEYS5:
        for p in allpoles[key]:
            for t in range(4):
                ys=[vec(k,key,p)[t] for k in range(7,12)]
                pred=lagr(list(range(7,12)),ys,kk)
                tot+=1; bad+= (pred!=vec(kk,key,p)[t])
    res[kk]=(tot,bad)
    print('k=',kk,'components',tot,'bad',bad)

print("\n== k=5 detail per key (components off the quartic)")
for key in KEYS5:
    off=[]
    for p in allpoles[key]:
        for t in range(4):
            ys=[vec(k,key,p)[t] for k in range(7,12)]
            pred=lagr(list(range(7,12)),ys,5)
            if pred!=vec(5,key,p)[t]: off.append((str(p[0])+'+'+str(p[1])+'r2',t,str(pred),str(vec(5,key,p)[t])))
    print(key,len(off),'off-quartic components at k=5')
    for o in off[:6]: print('   ',o)

print("\n== exact degree: 4th difference nonzero for some pole (k=7..11), per key; and 5-th zero")
for key in KEYS5:
    mx=0
    for p in allpoles[key]:
        seq=[vec(k,key,p) for k in range(7,12)]
        if any(c!=0 for c in fdiff(seq,4)): mx=4
    print(key,'degree-4 present:', mx==4)
