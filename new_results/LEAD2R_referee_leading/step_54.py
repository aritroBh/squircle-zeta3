import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 60
from math import factorial, comb
from lead2r_lib import e_coeffs_logE
E = e_coeffs_logE(16)
S = lambda u: mp.sqrt(1+u**2)
m_e = lambda u: (1-mp.sqrt(2)/S(u))/(1-u**2) if abs(1-u**2) > mp.mpf(10)**(-40) else mp.mpf(-1)/4
m_odd = lambda u: (1+mp.sqrt(2)*u/S(u))/(1-u**2)
m_diag = lambda u: (1-u)/2*m_e(u)
m_oddlast = lambda u: m_odd(u) - 1/(2*(1-u))
L1 = lambda u: mp.log((1+u)/2); L2 = lambda u: mp.log((1-u)/2)
def Mint(i, j, m):
    f = lambda u: L1(u)**i*L2(u)**j*m(u)
    tot = mp.mpf(0)
    # slices geometric toward u=-1 (x=(1+u)/2) and u=1
    x = mp.mpf(1)
    while x > mp.mpf(10)**(-50):
        tot += mp.quad(f, [-1+x/2, -1+x]); x /= 2
    y = mp.mpf(1)
    while y > mp.mpf(10)**(-50):
        tot += mp.quad(f, [1-y, 1-y/2]); y /= 2
    return tot
Mc = {}
def M(i,j,mname):
    m = {'e':m_e,'odd':m_odd,'diag':m_diag,'oddlast':m_oddlast}[mname]
    key = (i,j,mname)
    if key not in Mc: Mc[key] = Mint(i,j,m)
    return Mc[key]
worst = mp.mpf(0)
for d in range(1, 14):
    tot = mp.mpf(0)
    if d % 2 == 0:
        for j in range(d//2+1):
            mn = 'diag' if j == d//2 else 'e'
            tot += (-1)**(j+1)*comb(d,j)*M(d-j,j,mn)
    else:
        for j in range((d-1)//2+1):
            mn = 'oddlast' if j == (d-1)//2 else 'odd'
            tot += (-1)**(j+1)*comb(d,j)*M(d-j,j,mn)
    tgt = factorial(d)/2*E[d+1]
    r = abs(tot-tgt); worst = max(worst, r)
    print(f"d={d:2d} ({'even' if d%2==0 else 'odd'}): K_d - (d!/2)e_{{d+1}} = {mp.nstr(r,3)}")
print("worst:", mp.nstr(worst,3))
with open('progress.log','a') as f:
    f.write("step_21: full key sum eq(3) with 4 measures (incl diagonal/oddlast keys) == (d!/2)e_{d+1} for d=1..13 both parities, dps60: worst "+mp.nstr(worst,3)+"\n")
