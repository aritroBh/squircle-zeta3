import sys, pickle
sys.path.insert(0,'.')
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import mpmath as mp
mp.mp.dps = 50
import corner_series
from hpl import coeffs_mp

store = pickle.load(open('hpl_constructed_K7.pkl','rb'))
words = set(w for terms in store.values() for (e,w,q) in terms if len(w)>0)
print('distinct words:', len(words), flush=True)
N = 500
wc = {w: coeffs_mp(w, N) for w in words}
print('coeff arrays done', flush=True)
a = 2*mp.pi/3; c = mp.cos(a/2)
tot = mp.mpf(0)
for zkey, terms in store.items():
    zv = mp.mpf(1)
    for i,e in enumerate(zkey):
        if e: zv *= mp.zeta(i+2)**e
    s = mp.mpf(0)
    for (e,w,q) in terms:
        qq = mp.mpf(q.numerator)/q.denominator
        if len(w)==0:
            s += qq*c**e; continue
        ce = wc[w]; sv = mp.mpf(0); xp = c
        for m in range(1, N+1):
            sv += ce[m]*xp; xp *= c
        s += qq*c**e*sv
    tot += zv*s
Mc = corner_series.corner_coeffs(a, 7)
d = abs(tot - Mc[7])
print(f'k=7 HPL form (6553 terms) at c=1/2 vs corner_series M_7(2pi/3): residual {mp.nstr(d,3)} (~{mp.nstr(-mp.log10(d+mp.mpf("1e-999")),4)} digits)')
open('progress.log','a').write(f'step05e: k=7 HPL form == corner_series M_7(2pi/3): residual {mp.nstr(d,3)}\n')
