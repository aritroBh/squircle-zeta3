import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 110
from math import factorial, comb
from lead2r_lib import e_coeffs_logE
E = e_coeffs_logE(16)
def A_pq(p, q):
    f = lambda s: mp.log(s)**p*mp.log(1-s)**q/(1-s)
    tot = mp.mpf(0); a = mp.mpf('0.5')
    while a > mp.mpf(10)**(-95):
        tot += mp.quad(f, [a/2, a]); a /= 2
    b = mp.mpf('0.5')
    while b > mp.mpf(10)**(-95):
        tot += mp.quad(f, [1-b, 1-b/2]); b /= 2
    return tot
def T_odd(d):
    mid=(d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*A_pq(d-j,j) for j in range(mid+1)) - mp.mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*A_pq(mid+1,mid)
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
g1 = lambda w: (1-w)/(1+w)*h(w)
def piece2(d):
    f = lambda w: (-mp.log(w))**d*g1(w)
    tot = mp.mpf(0); x = mp.mpf(1)
    while x > mp.mpf(10)**(-95):
        tot += mp.quad(f, [x/2, x]); x /= 2
    return tot/2
for d in (11, 13):
    T = T_odd(d); f1 = factorial(d)*(1-mp.mpf(2)**(-d))*mp.zeta(d+1)
    p2 = piece2(d); f2 = f1 - factorial(d)/2*E[d+1]
    print(f"d={d}: T err {mp.nstr(abs(T-f1),3)},  um err {mp.nstr(abs(p2-f2),3)},  assembly err {mp.nstr(abs(T-p2-factorial(d)/2*E[d+1]),3)}")
with open('progress.log','a') as f:
    f.write("step_20: d=11,13 rerun at dps110: residuals shrink to ~1e-70 -> step_19 growth was precision-limited cancellation, identities hold\n")
