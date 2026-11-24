
import os, json
import mpmath as mp
mp.mp.dps = 60
base = "<project>"
h2 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n2_K80.json")))
eta2 = [mp.mpf(x) for x in h2['eta']]
d2 = [mp.mpf(x) for x in h2['d']]
# Just check: |a_n/n!|^{1/n} for the 2D tau series -> 2/pi (the 2D Borel radius is pi/2)
print("2D radius check |a_n/n!|^{1/n} -> 2/pi =", mp.nstr(2/mp.pi,12))
for n in [20,40,60,79]:
    v = abs(eta2[n]/mp.factorial(n))**(1/mp.mpf(n))
    print(f"  eta: n={n}: {mp.nstr(v,10)}")
# and check SRF1's n=2 step_37 finding A=2/pi: law form a_n/n! ~ C (2/pi)^n n^beta cos(n pi/2 + phi): fit beta
import math, numpy as np
b = [float(eta2[k]/mp.factorial(k)*(mp.pi/2)**k) for k in range(len(eta2))]
ks = list(range(30,80))
M = np.zeros((len(ks),4)); y=np.array([b[k] for k in ks])
for i,k in enumerate(ks):
    M[i,0]=k**-1.5*math.cos(math.pi/2*k); M[i,1]=k**-1.5*math.sin(math.pi/2*k)
    M[i,2]=k**-2.5*math.cos(math.pi/2*k); M[i,3]=k**-2.5*math.sin(math.pi/2*k)
x,*_ = np.linalg.lstsq(M,y,rcond=None)
rms = float(np.sqrt(np.mean((y-M@x)**2)))
print("2D fit beta=-1.5: C=",math.hypot(x[0],x[1]),"phi=",math.atan2(x[1],x[0]),"rms=",rms)
# paper law: f_n/n! ~ sqrt2 cos(pi(3-2n)/4) (2/pi)^n n^-3/2 -> b_k should be sqrt2 cos(pi(3-2k)/4)
for k in [40,60,79]:
    print("  b_k=",b[k], " paper law:", mp.sqrt(2)*mp.cos(mp.pi*(3-2*k)/4))
