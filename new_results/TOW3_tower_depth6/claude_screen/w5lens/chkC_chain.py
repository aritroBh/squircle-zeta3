# chkC_chain.py -- end-to-end chain from TOW3's PROVED closed forms (not data) to the degree-4 kernels:
#  closed-form towers u(delta)(k) = e00 * sum_s C(k-5,s) e(delta,s)(d=5,m)  for keys (5,0),(4,1),(3,2)
#  -> (a) compare with exact TOWN1 towers (k=5..52, delta<=6)
#  -> (b) feed into the Lemma-2 map (eq:triangular at degree 4: W40=-5G50 w1-G41 w2, W31=-4G41 w1-2G32 w2, W22=-3G32 w1)
#  -> (c) compare with exact reducer kernels k=5..21; (d) cubic law, no +-i simple pole, zero const for k=5..60.
import os, sys, pickle, time
from fractions import Fraction as F
sys.dont_write_bytecode=True
R='<project>/research/08_theory/'
sys.path.insert(0, R+'WTS_weight_k_structure/theory')
sys.path.insert(0, R+'TOWN1_tower_large_n')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import tow3lib as T3
import kernel_from_tower as kft
from fastred import QSI, ZERO, pf_add, pf_scale, pf_mul_DL, DL1, DL2
from math import comb

H=pickle.load(open(R+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
KEYS=[(5,0),(4,1),(3,2)]
# relative e(delta,s) at (d,m)=(5,m) as Gaussian rationals (x,y)
def gauss(expr):
    expr=sp.expand(expr); re_,im_=expr.as_real_imag()
    re_=sp.nsimplify(re_) if not re_.is_Rational else re_
    im_=sp.nsimplify(im_) if not im_.is_Rational else im_
    assert re_.is_Rational and im_.is_Rational, expr
    return F(int(re_.p),int(re_.q)), F(int(im_.p),int(im_.q))
E={}
for key in KEYS:
    m=key[0]
    for delta in range(0,7):
        for s in range(0,delta//2+1):
            E[(key,delta,s)]=gauss(T3.e(delta,s,5,m))
print('relative e(delta,s) evaluated at d=5 for the three keys, count',len(E))
e00={key:QSI(*H[key][0][5]) for key in KEYS}
def u_closed(key,delta,k):
    x=F(0);y=F(0)
    for s in range(0,delta//2+1):
        c=comb(k-5,s) if k>=5 else 0
        ex,ey=E[(key,delta,s)]
        x+=c*ex; y+=c*ey
    return e00[key]*QSI(x,0,y,0)
# (a) closed-form towers vs exact towers
bad=0;n=0
for key in KEYS:
    for delta in range(0,7):
        for k in range(5,53):
            n+=1
            if tuple(u_closed(key,delta,k).__getattribute__(a) for a in ('p0','p1','p2','p3'))!=tuple(H[key][delta][k]):
                bad+=1
                if bad<6: print('tower mismatch',key,delta,k)
print(f'(a) closed-form degree-5 towers vs exact TOWN1 towers: {n} cells (3 keys x depths 0..6 x k=5..52), {bad} mismatches')
# (b) Lemma-2 map from CLOSED-FORM towers
def h_orders_closed(key,k):
    return {8-delta:u_closed(key,delta,k) for delta in range(0,7)}
def kernel_deg4_from(hfun,k):
    g50,g41,g32=(kft.gm_from_tower(hfun(key,k)) for key in KEYS)
    W40=pf_add(pf_scale(pf_mul_DL(g50,DL1),F(-5)),pf_scale(pf_mul_DL(g41,DL2),F(-1)))
    W31=pf_add(pf_scale(pf_mul_DL(g41,DL1),F(-4)),pf_scale(pf_mul_DL(g32,DL2),F(-2)))
    W22=pf_scale(pf_mul_DL(g32,DL1),F(-3))
    out={}
    for key,W in [((4,0),W40),((3,1),W31),((2,2),W22)]:
        simple,higher,const=kft.pf_split(W)
        out[key]=(simple,const)
    return out
def flat(simple,const):
    """flatten simple part: real residues (4 comps), ipos/ineg all entries, const"""
    d={}
    for r,c in simple[0].items():
        for nm,v in zip(('p0','p1','p2','p3'),(c.p0,c.p1,c.p2,c.p3)):
            if v!=0: d[('real',(r.p0,r.p1),nm)]=v
    for tag,idx in (('ipos',1),('ineg',2)):
        for j,c in enumerate(simple[idx]):
            for nm,v in zip(('p0','p1','p2','p3'),(c.p0,c.p1,c.p2,c.p3)):
                if v!=0: d[(tag,j,nm)]=v
    for nm,v in zip(('p0','p1','p2','p3'),(const.p0,const.p1,const.p2,const.p3)):
        if v!=0: d[('const',nm)]=v
    return d
KC={}
t0=time.time()
for k in range(5,61):
    KC[k]={key:flat(*v) for key,v in kernel_deg4_from(h_orders_closed,k).items()}
print('kernels from closed-form towers built for k=5..60 in %.1fs'%(time.time()-t0))
# (c) vs exact reducer kernels k=5..21
bad=0;n=0
for k in range(5,22):
    fn=R+f'WTS_weight_k_structure/theory/kern_{k}.pkl' if k<=18 else R+f'TWOEND1_weight_poles/red/kern_{k}.pkl'
    Kref=pickle.load(open(fn,'rb'))
    for key in [(4,0),(3,1),(2,2)]:
        ref={}
        P=Kref[key]
        for pole,co in P['real'].items():
            for nm,v in zip(('p0','p1','p2','p3'),co):
                if v!=0: ref[('real',pole,nm)]=v
        for tag in ('ipos','ineg'):
            for j,co in enumerate(P[tag]):
                for nm,v in zip(('p0','p1','p2','p3'),co):
                    if v!=0: ref[(tag,j,nm)]=v
        for nm,v in zip(('p0','p1','p2','p3'),P['const']):
            if v!=0: ref[('const',nm)]=v
        n+=1
        if ref!=KC[k][key]:
            bad+=1; print('kernel mismatch',k,key)
print(f'(c) kernels from CLOSED-FORM towers vs exact reducer kernels: {n} (k=5..21 x 3 keys), {bad} mismatches')
# (d) cubic law and absence of +-i / const, k=5..60
bad_i=0
for k in range(5,61):
    for key in [(4,0),(3,1),(2,2)]:
        for c in KC[k][key]:
            if c[0] in ('ipos','ineg','const'): bad_i+=1
            if c[0]=='real' and c[2] in ('p2','p3'): bad_i+=1
print('(d1) +-i simple-pole / const / imaginary-residue entries over k=5..60 (should be 0):',bad_i)
nz=0;tot=0;deg3=0
for key in [(4,0),(3,1),(2,2)]:
    cs=set().union(*[set(KC[k][key]) for k in range(5,61)])
    for c in cs:
        seq=[KC[k][key].get(c,F(0)) for k in range(5,61)]
        tot+=1
        d4=seq
        for _ in range(4): d4=[d4[i+1]-d4[i] for i in range(len(d4)-1)]
        if any(x!=0 for x in d4): nz+=1
        d3=seq
        for _ in range(3): d3=[d3[i+1]-d3[i] for i in range(len(d3)-1)]
        if any(x!=0 for x in d3): deg3+=1
print(f'(d2) components (key,pole,QSI-part) checked {tot}; 4th difference nonzero (k=5..60): {nz}; components of exact degree 3: {deg3}')
