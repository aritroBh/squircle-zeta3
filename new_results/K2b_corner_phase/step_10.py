
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
for key in ['pi_3','pi_2','2pi_3','3pi_4','5pi_6']:
    a = ALPHA[key]; d = load(key)
    A = mp.sqrt(2*a*mp.sin(a))
    print('===', key, 'A=', mp.nstr(A,12))
    for k in [50,100,150,200,250,300]:
        N = k+1
        th = mp.sqrt(2*a*N)
        W = -d[k]*N*mp.e**th/A   # expect ~ Re[e^{i th}(1 + k1/sqrtN + ...)] = cos th + ...
        print('  k=%3d  W=%s  cos(th)=%s' % (k, mp.nstr(mp.re(W),16), mp.nstr(mp.cos(th),10)))
