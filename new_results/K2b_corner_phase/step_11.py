
wd = 'research/08_theory/K2b_corner_phase/'
s = open(wd+'k2utils.py').read()
s = s.replace("return mpmath.mpf(s[:idx]), mpmath.mpf(s[idx:].replace('j','').replace('i',''))",
              "return mpmath.mpf(s[:idx].replace(' ','')), mpmath.mpf(s[idx:].replace('j','').replace('i','').replace(' ',''))")
open(wd+'k2utils.py','w').write(s)
import sys; sys.path.insert(0,wd)
import k2utils; importlib_ok=True
from k2utils import *
mp.dps = 80
for key in ['pi_3','pi_2','2pi_3','3pi_4','5pi_6']:
    a = ALPHA[key]; d = load(key)
    A = mp.sqrt(2*a*mp.sin(a))
    print('===', key, 'A=', mp.nstr(A,12))
    for k in [50,100,150,200,250,300]:
        N = k+1
        th = mp.sqrt(2*a*N)
        W = -d[k]*N*mp.e**th/A
        print('  k=%3d  W=%s  cos(th)=%s' % (k, mp.nstr(mp.re(W),16), mp.nstr(mp.cos(th),10)))
