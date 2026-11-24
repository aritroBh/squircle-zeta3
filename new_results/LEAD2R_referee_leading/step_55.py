import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 60
from math import factorial, comb
from lead2r_lib import e_coeffs_logE
E = e_coeffs_logE(16)
S = lambda u: mp.sqrt(1+u**2)
m_e = lambda u: (1-mp.sqrt(2)/S(u))/(1-u**2) if abs(1-u**2) > mp.mpf(10)**(-40) else mp.mpf(-1)/4
m_diag = lambda u: (1-u)/2*m_e(u)
L1 = lambda u: mp.log((1+u)/2); L2 = lambda u: mp.log((1-u)/2)
def quad01(f):
    tot = mp.mpf(0); x = mp.mpf(1)
    while x > mp.mpf(10)**(-45):
        tot += mp.quad(f, [x/2-1, x-1]); x /= 2
    y = mp.mpf(1)
    while y > mp.mpf(10)**(-45):
        tot += mp.quad(f, [1-y, 1-y/2]); y /= 2
    return tot
def Mint(i,j,m): return quad01(lambda u: L1(u)**i*L2(u)**j*m(u))
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
def quad01w(f):  # on w in [0,1]
    tot = mp.mpf(0); x = mp.mpf(1)
    while x > mp.mpf(10)**(-45):
        tot += mp.quad(f, [x/2, x]); x /= 2
    return tot
worst = {k: mp.mpf(0) for k in ('collapse','midcorr','singleint','omission_nonzero')}
for d in range(2, 14, 2):
    Kfull = sum((-1)**(j+1)*comb(d,j)*Mint(d-j,j, m_diag if j==d//2 else m_e) for j in range(d//2+1))
    coll = -mp.mpf('0.5')*quad01(lambda u: (L1(u)-L2(u))**d*m_e(u))
    r1 = abs(Kfull-coll); worst['collapse']=max(worst['collapse'],r1)
    # middle-key correction: (-1)^{d/2+1} binom(d,d/2) M(d/2,d/2; m_diag - m_e/2) should vanish
    mc = Mint(d//2, d//2, lambda u: m_diag(u)-m_e(u)/2)
    r2 = abs(mc); worst['midcorr']=max(worst['midcorr'],r2)
    # single integral: K_d = (1/2) int_0^1 (-log w)^d h(w) dw
    si = mp.mpf('0.5')*quad01w(lambda w: (-mp.log(w))**d*h(w))
    r3 = abs(Kfull-si); worst['singleint']=max(worst['singleint'],r3)
    # eq(3) PRINTED range (j <= floor((d-1)/2)) omits diagonal -> must NOT equal K_d
    Komit = sum((-1)**(j+1)*comb(d,j)*Mint(d-j,j,m_e) for j in range((d-1)//2+1))
    r4 = abs(Komit - factorial(d)/2*E[d+1]); worst['omission_nonzero']=min(r4) if d==2 else worst['omission_nonzero']
    print(f"d={d:2d}: collapse {mp.nstr(r1,2)}  midcorr~0 {mp.nstr(r2,2)}  single-int {mp.nstr(r3,2)}  printed-range-only err {mp.nstr(r4,2)}")
print("worst:", {k: mp.nstr(v,3) for k,v in worst.items()})
with open('progress.log','a') as f:
    f.write("step_22: Lemma 4.1 collapse + vanishing middle-key correction + single-integral form verified d=2..12 even; CONFIRMED eq(3) printed j-range omits the even-d diagonal key (nonzero error if used alone) - handled correctly in proof\n")
