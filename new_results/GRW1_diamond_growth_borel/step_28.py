
import mpmath as mp
mp.mp.dps = 110
import json, numpy as np

P = json.load(open("pade9090.json"))
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
den_d = np.array([complex(v) for v in den[::-1]])
rts = np.roots(den_d)
cand = sorted(rts, key=abs)[:24]
def denf(z): return sum(den[k]*z**k for k in range(len(den)))
refined = []
for z0 in cand:
    try:
        z = mp.findroot(denf, mp.mpc(z0.real, z0.imag), tol=mp.mpf('1e-60'), maxsteps=60)
        refined.append(z)
    except Exception as e:
        pass
refined.sort(key=abs)
print("refined poles of Pade-Borel [90/90] (smallest |z|); pi/2 =", mp.nstr(mp.pi/2, 12))
for z in refined[:12]:
    print(f"   {mp.nstr(z, 25)}   |z|={mp.nstr(abs(z), 15)}  |z|/(pi/2)-1={mp.nstr(abs(z)/(mp.pi/2)-1, 4)}")
print("any pole with |z| < pi/2:", any(abs(z) < mp.pi/2 for z in refined))
print("any pole within 0.05 of the real axis (|Im|<0.05, |z|<20):", any(abs(mp.im(z))<0.05 and abs(z)<20 for z in refined))
