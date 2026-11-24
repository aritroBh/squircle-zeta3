
import sys; sys.path.insert(0,'research/08_theory/K2b_corner_phase/')
from k2utils import *
mp.dps = 80
for key in ['pi_3','pi_2','2pi_3','3pi_4','5pi_6']:
    a = ALPHA[key]; d = load(key)
    A = mp.sqrt(2*a*mp.sin(a))
    print('===', key, 'alpha/pi=', mp.nstr(a/mp.pi,5), 'A=', mp.nstr(A,12))
    for k in [50,100,150,200,250,300]:
        N = k+1
        th = mp.sqrt(2*a*N)
        W = -d[k]*N*mp.e**th/A   # should ~ cos(th) + Re[e^{i th} k1]/sqrt(N)+...
        print('  k=%3d  W=%s  cos(th)=%s' % (k, mp.nstr(mp.re(W),18), mp.nstr(mp.cos(th),12)))
