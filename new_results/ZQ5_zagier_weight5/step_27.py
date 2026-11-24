import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib as W
P = 830
mp.mp.dps = P+25
x = mp.mpf(1)/2
# manual with counting
s = mp.mpf(0); k = 1; xk = x
thresh = mp.mpf(10)**(-(P+15))
while True:
    t = xk/mp.power(k,5)
    s += t
    if abs(t) < thresh and k > 30: break
    xk *= x; k += 1
print("manual terms:", k, "last |t|:", mp.nstr(abs(t),3))
v = W.li5_series(x, P)
print("polylog - manual:", mp.nstr(mp.polylog(5,x) - s, 3))
print("polylog - li5_series:", mp.nstr(mp.polylog(5,x) - v, 3))
# check thresh
print("thresh:", mp.nstr(thresh, 3))
# when does 2^-k/k^5 first drop below 1e-845?
import math
kk = 1
xk2 = mp.mpf('0.5')
while xk2/mp.power(kk,5) > thresh:
    xk2 *= x; kk += 1
print("first k below thresh:", kk)
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',src],capture_output=True,text=True,timeout=600)
print(r.stdout); print("STDERR:", r.stderr[-800:])
