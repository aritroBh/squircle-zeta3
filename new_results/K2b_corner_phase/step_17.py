
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
d = load('pi_2'); a = mp.pi/2; A = mp.sqrt(mp.pi)
# D13b convention: M = -(sqrt(pi)/N) e^-th [ cos th (1 + cr/sqrtN) + (ci/sqrtN) sin th ]
cr = 1/(2*mp.sqrt(mp.pi)) + mp.sqrt(mp.pi)/4 + mp.pi**mp.mpf('1.5')/48
ci = -1/(2*mp.sqrt(mp.pi)) + mp.sqrt(mp.pi)/4 - mp.pi**mp.mpf('1.5')/48
print('cr,ci =', mp.nstr(cr,18), mp.nstr(ci,18))
print('k   rho_D13b*N   (should be bounded/decay)   rho_fit*N with k1=0.7776')
for k in [60,80,100,150,200,250,300]:
    N = mp.mpf(k+1); th = mp.sqrt(mp.pi*N)
    W = -d[k]*N*mp.e**th/A
    r1 = W - (mp.cos(th)*(1+cr/mp.sqrt(N)) + ci/mp.sqrt(N)*mp.sin(th))
    r2 = W - (mp.cos(th)*(1+mp.mpf('0.777555688923117')/mp.sqrt(N)))
    print(k, mp.nstr(r1*N,10), '   ', mp.nstr(r2*N,10))
