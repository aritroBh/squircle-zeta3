# chkB_kernels.py -- own checks (2),(3),(4) from the exact degree-4 kernels kern_k.pkl (k<=18) + kern_19..21.pkl
import pickle, sys, json
from fractions import Fraction as Fr
R='<project>/research/08_theory/'
sys.path.insert(0, R+'WTS_weight_k_structure/theory')
KEYS=[(4,0),(3,1),(2,2)]
K={}
for k in range(5,22):
    fn = R+f'WTS_weight_k_structure/theory/kern_{k}.pkl' if k<=18 else R+f'TWOEND1_weight_poles/red/kern_{k}.pkl'
    K[k]=pickle.load(open(fn,'rb'))
Z4=(Fr(0),)*4
POLES_ALLOWED={(Fr(0),Fr(0)),(Fr(1),Fr(1)),(Fr(-1),Fr(-1)),(Fr(1),Fr(-1))}   # 0, a, -a, -b
names={(Fr(0),Fr(0)):'0',(Fr(1),Fr(1)):'a',(Fr(-1),Fr(-1)):'-a',(Fr(1),Fr(-1)):'-b'}
def comp(k,key):
    """flatten kernel to dict {('real',pole,c):value, ('ipos',idx,c):..., ('const',c):...}"""
    P=K[k][key]; out={}
    for pole,co in P['real'].items():
        assert pole in POLES_ALLOWED, ('unexpected pole',k,key,pole)
        for c in range(4): out[('real',pole,c)]=co[c]
    for idx,co in enumerate(P['ipos']):
        for c in range(4): out[('ipos',idx,c)]=co[c]
    for idx,co in enumerate(P['ineg']):
        for c in range(4): out[('ineg',idx,c)]=co[c]
    for c in range(4): out[('const',c)]=P['const'][c]
    return out
# ---------- (2) no +-i simple pole, const zero, only real poles in {0,a,-a,-b}, c2=c3=0 on real residues
n2=0; bad2=0
for k in range(5,22):
    for key in KEYS:
        P=K[k][key]
        n2+=1
        if P['ipos'] or P['ineg'] or any(x!=0 for x in P['const']): bad2+=1; print('(2) BAD',k,key,P['ipos'],P['ineg'],P['const'])
        for pole,co in P['real'].items():
            if co[2]!=0 or co[3]!=0: bad2+=1; print('(2) imaginary part in real residue',k,key,pole,co)
print(f'(2) kernel-levels checked {n2} (k=5..21 x 3 keys): ipos/ineg lists nonempty or const nonzero or imag part: {bad2} bad')
# also list the sets of poles present
for key in KEYS:
    ps=set()
    for k in range(5,22):
        ps|={names[p] for p,co in K[k][key]['real'].items() if any(co)}
    print('   key',key,'pole set over k=5..21:',sorted(ps))
# ---------- (3) K5 = 4K6-6K7+4K8-K9 and cubic through 6..9 reproduces 10..21
allk=lambda key,ks:set().union(*[set(comp(k,key)) for k in ks])
nc=0; bad3=0
for key in KEYS:
    cs=allk(key,range(5,22))
    C={k:comp(k,key) for k in range(5,22)}
    for c in sorted(cs,key=str):
        g=lambda k: C[k].get(c,Fr(0))
        nc+=1
        if g(5)!=4*g(6)-6*g(7)+4*g(8)-g(9): bad3+=1; print('(3) K5 relation FAILS',key,c)
        # Lagrange cubic through 6..9
        for k in range(10,22):
            val=Fr(0)
            for a in range(6,10):
                li=Fr(1)
                for b in range(6,10):
                    if a!=b: li*=Fr(k-b,a-b)
                val+=g(a)*li
            if val!=g(k): bad3+=1; print('(3) cubic(6..9) fails at k=',k,key,c)
        # 4th differences all vanishing on k=5..21 (13 windows)
        for k in range(5,18):
            fd=sum((-1)**j*__import__('math').comb(4,j)*g(k+j) for j in range(5))
            if fd!=0: bad3+=1; print('(3) 4th diff nonzero',k,key,c)
print(f'(3) components tested (key,pole/imag/const,QSI comp): {nc}; failures: {bad3}')
# exact degree 3: the third difference is nonzero for at least one component per key; k^3 coefficient
for key in KEYS:
    cs=allk(key,range(5,22)); C={k:comp(k,key) for k in range(5,22)}
    nz=[c for c in cs if (C[8].get(c,0)-3*C[7].get(c,0)+3*C[6].get(c,0)-C[5].get(c,0))!=0]
    print('   key',key,'#components with nonzero 3rd difference:',len(nz))
# ---------- (4) differences vs WT5 data
D=json.load(open(R+'WT5_weight5_all_k/deltas_deg4.json'))
KD=json.load(open(R+'WT5_weight5_all_k/kernels_deg4.json'))
def fr(s): return Fr(s)
inv={'0':(Fr(0),Fr(0)),'a':(Fr(1),Fr(1)),'-a':(Fr(-1),Fr(-1)),'-b':(Fr(1),Fr(-1))}
def from_json(lst):
    out={}
    for nm,(c0,c1) in lst:
        out[('real',inv[nm],0)]=fr(c0); out[('real',inv[nm],1)]=fr(c1)
    return out
def diffk(key,j,base=5):
    C={k:comp(k,key) for k in range(base,base+j+1)}
    cs=set().union(*[set(C[k]) for k in C])
    out={}
    for c in cs:
        out[c]=sum((-1)**(j-t)*__import__('math').comb(j,t)*C[base+t].get(c,Fr(0)) for t in range(j+1))
    return out
bad4=0; n4=0
for ks in ('4,0','3,1','2,2'):
    key=tuple(int(x) for x in ks.split(','))
    for j,nm in ((1,'D1'),(2,'D2'),(3,'D3')):
        mine=diffk(key,j)
        theirs=from_json(D[ks][nm])
        cs=set(mine)|set(theirs)
        for c in cs:
            n4+=1
            if mine.get(c,Fr(0))!=theirs.get(c,Fr(0)): bad4+=1; print('(4) mismatch',key,nm,c,mine.get(c),theirs.get(c))
print(f'(4) difference kernels at k=5 vs WT5 deltas_deg4.json: {n4} components compared, {bad4} mismatches')
# also kernels_deg4.json k=5..14 vs pickles
bad4b=0; n4b=0
for ksn,kd in KD.items():
    k=int(ksn)
    for ks,lst in kd.items():
        key=tuple(int(x) for x in ks.split(','))
        th=from_json(lst); mine=comp(k,key)
        for c in set(th)|{c for c in mine if c[0]=='real'}:
            n4b+=1
            if mine.get(c,Fr(0))!=th.get(c,Fr(0)): bad4b+=1; print('(4b) mismatch',k,key,c)
print(f'(4b) kernels_deg4.json (k=5..14) vs kern pickles: {n4b} components, {bad4b} mismatches')
# print the actual K5..K9 for the record
for key in KEYS:
    for k in (5,6):
        P=K[k][key]['real']
        print(key,k,{names[p]:tuple(str(x) for x in co[:2]) for p,co in sorted(P.items(),key=lambda t:names[t[0]])})
