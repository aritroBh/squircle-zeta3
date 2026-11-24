# chkE_rho.py -- initial data rho(delta)=a_delta*rho(0) at k=d, delta=0..6, all keys of degree d=5..16 (incl. q=0 and diagonal), vs exact towers
import sys, os, pickle
from fractions import Fraction as F
sys.dont_write_bytecode=True
R='<project>/research/08_theory/'
sys.path.insert(0, R+'WTS_weight_k_structure/theory'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp, tow3lib as T3
from fastred import QSI
H=pickle.load(open(R+'TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')
def gauss(expr):
    expr=sp.expand(expr); re_,im_=expr.as_real_imag()
    return F(int(sp.Rational(re_).p),int(sp.Rational(re_).q)),F(int(sp.Rational(im_).p),int(sp.Rational(im_).q))
n=bad=0
for d in range(5,17):
    for p_ in range((d+1)//2, d+1):
        q_=d-p_
        key=(p_,q_)
        r0=QSI(*H[key][0][d])
        for delta in range(1,7):
            if delta>2*d-4: continue
            x,y=gauss(T3.arho(delta).subs({T3.d:d,T3.m:p_}))
            pred=r0*QSI(x,0,y,0)
            n+=1
            if (pred.p0,pred.p1,pred.p2,pred.p3)!=tuple(H[key][delta][d]):
                bad+=1; print('rho mismatch',key,delta)
print(f'rho(delta)=a_delta*rho(0) at k=d, delta=1..6, keys of degree 5..16 (all p>=q incl q=0 and diagonal): {n} cells, {bad} mismatches')
