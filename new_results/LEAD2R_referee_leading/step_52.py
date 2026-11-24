import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 55
from math import factorial, comb
from lead2r_lib import e_coeffs_logE
E = e_coeffs_logE(16)
def A_pq(p, q):
    f = lambda s: mp.log(s)**p*mp.log(1-s)**q/(1-s)
    tot = mp.mpf(0); a = mp.mpf('0.5')
    while a > mp.mpf(10)**(-40):
        tot += mp.quad(f, [a/2, a]); a /= 2
    b = mp.mpf('0.5')
    while b > mp.mpf(10)**(-40):
        tot += mp.quad(f, [1-b, 1-b/2]); b /= 2
    return tot
Acache = {}
def A(p,q):
    if (p,q) not in Acache: Acache[(p,q)] = A_pq(p,q)
    return Acache[(p,q)]
def T_odd(d):
    mid=(d-1)//2
    return sum((-1)**(j+1)*comb(d,j)*A(d-j,j) for j in range(mid+1)) - mp.mpf('0.5')*(-1)**(mid+1)*comb(d,mid)*A(mid+1,mid)
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
g1 = lambda w: (1-w)/(1+w)*h(w)
def piece2(d):  # (u m_e part) = (1/2) int_0^1 (-log w)^d g1(w) dw  (odd d; from the proved reduction)
    f = lambda w: (-mp.log(w))**d*g1(w)
    tot = mp.mpf(0); x = mp.mpf(1)
    while x > mp.mpf(10)**(-40):
        tot += mp.quad(f, [x/2, x]); x /= 2
    return tot/2
print(f"{'d':>3} {'T_d/d! - (1-2^-d)zeta(d+1)':>28} {'umpart/d! - rhs':>28} {'(T-umpart)/d! - e/2':>28}")
worst = {k: mp.mpf(0) for k in ('T','um','asm')}
for d in range(1, 14, 2):
    T = T_odd(d); f1 = factorial(d)*(1-mp.mpf(2)**(-d))*mp.zeta(d+1)
    p2 = piece2(d); f2 = f1 - factorial(d)/2*E[d+1]
    asm = T - p2; tgt = factorial(d)/2*E[d+1]
    r1, r2, r3 = abs(T-f1), abs(p2-f2), abs(asm-tgt)
    worst['T']=max(worst['T'],r1); worst['um']=max(worst['um'],r2); worst['asm']=max(worst['asm'],r3)
    print(f"{d:3d} {mp.nstr(r1,3):>28} {mp.nstr(r2,3):>28} {mp.nstr(r3,3):>28}")
print("worst:", {k: mp.nstr(v,3) for k,v in worst.items()})
with open('progress.log','a') as f:
    f.write("step_19: Thm 5.4 T_d=d!(1-2^-d)zeta(d+1) + Lemma 5.1 (7) + assembly (Thm 5.5) for odd d=1..13 via quadrature dps55: worst T "+mp.nstr(worst['T'],3)+", um "+mp.nstr(worst['um'],3)+", assembly "+mp.nstr(worst['asm'],3)+"\n")
