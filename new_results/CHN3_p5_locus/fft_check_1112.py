"""Cauchy/FFT check of c11, c12 at a point, using h(t) representation:
Vol/32 = u3 int_0^1 t^(u3-1) h(t) dt, h(t)=(1-t)^(u2+u4) 2F1(-u1,u2;1+u2;1-t) 2F1(-u5,u4;1+u4;1-t),
analytically continued as h(0) + u3 int t^(u3-1)(h(t)-h(0)) dt for small complex u.
Compare with engine prediction (point values from forms1112 or scalar engine)."""
import sys, os, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
import mpmath as mp
import mzv_vals
POINT=(1,2,3,4,5)
def vol(a,u):
    u1,u2,u3,u4,u5=[x*u for x in a]
    h=lambda t:(1-t)**(u2+u4)*mp.hyp2f1(-u1,u2,1+u2,1-t)*mp.hyp2f1(-u5,u4,1+u4,1-t)
    h0=mp.hyp2f1(-u1,u2,1+u2,1)*mp.hyp2f1(-u5,u4,1+u4,1)
    return h0+u3*mp.quad(lambda t:t**(u3-1)*(h(t)-h0),[0,mp.mpf('0.001'),mp.mpf('0.05'),mp.mpf('0.3'),1])
def cauchy(a, W, N=96, dps=70):
    R=1/mp.mpf(max(a[0],a[4],a[1]+a[3],a[0]+a[2],a[2]+a[4])); r=R/6
    mp.mp.dps=dps
    vals={}
    for j in range(N//2+1):
        vals[j]=vol(a, r*mp.expjpi(mp.mpf(2*j)/N))
    for j in range(N//2+1,N): vals[j]=mp.conj(vals[N-j])
    return mp.fsum(vals[j]*mp.expjpi(-mp.mpf(2*W*j)/N) for j in range(N))/N/r**W
if __name__=='__main__':
    import mzv_exact as M
    for W in (11,12):
        cW = cauchy(POINT, W)
        # engine prediction from test_pt_w1112.pkl
        d = pickle.load(open('test_pt_w1112.pkl','rb'))
        pred = d['c11'] if W==11 else d['c12']
        T = M.DSTableFlint(W) if W>=11 else M.DSTable(W)
        # build numerical basis from monomials
        zv = {2:mp.zeta(2),3:mp.zeta(3),5:mp.zeta(5),7:mp.zeta(7),9:mp.zeta(9),11:mp.zeta(11)}
        import itertools
        def mzv_of(word):
            if len(word)==1: return mp.zeta(word[0])
            if len(word)==2: return mzv_vals.dz(word[0],word[1])
            raise ValueError(word)
        def monoval(mon):
            v=mp.mpf(1)
            for L in mon: v*=mzv_of(L)
            return v
        basis_vals=[]
        for mon in T.monomials:
            # mon is tuple of words
            basis_vals.append(monoval(mon))
        pv = mp.fsum(mp.mpf(c.numerator)/c.denominator*b for c,b in zip(pred,basis_vals))
        print(f"W={W}: cauchy={mp.nstr(cW.real,30)} imag={mp.nstr(cW.imag,2)}")
        print(f"   engine={mp.nstr(pv,30)}")
        print(f"   diff={mp.nstr(cW.real-pv,3)}")
