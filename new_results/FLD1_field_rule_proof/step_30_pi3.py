"""step30_pi3.py -- high-precision k=4 at pi/3: HPL form vs corner_series, and vs G2's B4 value."""
import sys, pickle
import mpmath as mp
sys.path.insert(0,'.')
from hpl import coeffs_mp
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import corner_series
mp.mp.dps = 110
store = pickle.load(open('hpl_constructed_K5.pkl','rb'))
def zeta_of(zkey):
    v=mp.mpf(1)
    for i,e in enumerate(zkey):
        if e: v*=mp.zeta(i+2)**e
    return v
N=2500
words=set()
for (k,zk),terms in store.items():
    if k!=4: continue
    for (e,w,q) in terms:
        if len(w)>0: words.add(w)
wc={w: coeffs_mp(w,N) for w in words}
def Mk4(c):
    tot=mp.mpf(0)
    for (kk,zk),terms in store.items():
        if kk!=4: continue
        s=mp.mpf(0)
        for (e,w,q) in terms:
            qq=mp.mpf(q.numerator)/q.denominator
            if len(w)==0: s+=qq*c**e; continue
            ce=wc[w]; sv=mp.mpf(0); xp=c
            for m in range(1,N+1): sv+=ce[m]*xp; xp*=c
            s+=qq*c**e*sv
        tot+=zeta_of(zk)*s
    return tot
a=mp.pi/3; c=mp.cos(a/2)
v = Mk4(c)
Mc = corner_series.corner_coeffs(a,4)[4]
d = abs(v-Mc)
print(f'pi/3 (c={mp.nstr(c,10)}):')
print('  HPL form M4  =', mp.nstr(v, 60))
print('  corner_series=', mp.nstr(Mc, 60))
print(f'  |diff| = {mp.nstr(d,3)}  (~{mp.nstr(-mp.log10(d+mp.mpf("1e-999")),4)} digits)')
# Taylor-convergence: change N
