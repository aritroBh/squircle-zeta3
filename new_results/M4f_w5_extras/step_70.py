T='<project>/research/08_theory'
txt = open(f'{T}/M4a_decomposition/M4_800.txt').read()
print(repr(txt[:400]))
print('...')
import re
val = re.sub(r'\s','',txt).split('=')[-1]
print(repr(val[:80]), 'len', len(val))
# M4c's W34 formula at 880 dps
import mpmath as mp
mp.mp.dps = 880
r2, lam, L = mp.sqrt(2), mp.log(2), mp.log(1+mp.sqrt(2))
b = r2-1
W3p = (7*r2/32 - mp.mpf(3)/16)*mp.zeta(3) + (1+r2)*mp.polylog(3,b) + (1-r2)*mp.polylog(3,-b) \
      - 2*r2*mp.polylog(3,1/r2) + (2*L - r2*lam)*mp.polylog(2,b) \
      + (r2/12 - mp.mpf(1)/8)*mp.pi**2*L - (1+r2)/24*mp.pi**2*lam + (1+r2)/24*lam**3 \
      + r2/4*lam**2*L - (1+r2)/2*lam*L**2 + (5-r2)/6*L**3
W34 = W3p - mp.zeta(3)/8 - mp.pi**2*lam/24 - lam**3/12
print('M4c W34 =', mp.nstr(W34, 30))
