import os, json, time
os.chdir("research/08_theory/K4b_cubo_finish")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
import importlib, k4_engine, k4_quad
importlib.reload(k4_engine); importlib.reload(k4_quad)
from k4_engine import *
from k4_quad import *

WMIN = mpf('1e-68'); LMIN = mpf('1e-60')
MT = fast_tetra(K=4)
KK = 4

def regionB_vec(ph, HI, TI):
    thB = th_B(ph)
    cph, sph = mp.cos(ph), mp.sin(ph)
    sig = SIG_B
    def g(th):
        sth, cth = mp.sin(th), mp.cos(th)
        ux, uy, uz = sth*cph, sth*sph, cth
        for j in range(4):
            aj = sig[j]*(MT.n[j][0]*ux + MT.n[j][1]*uy + MT.n[j][2]*uz)
            if aj <= 0:
                return None
        Fs = MT.F_series(ux, uy, uz, sig)
        return [Fs[k]*sth for k in range(KK+1)]
    acc = [mpf(0)]*(KK+1)
    for (x,w) in tanh_sinh_nodes(mpf(0), thB, HI, TI, WMIN):
        gx = g(x)
        if gx is None: continue
        for k in range(KK+1): acc[k] += w*gx[k]
    return acc

t0=time.time()
S_B128 = [mpf(0)]*(KK+1)
for (x,w) in tanh_sinh_nodes(mp.pi/4, mp.pi/2, mpf(1)/32, mpf('4.7'), WMIN):
    vB = regionB_vec(x, mpf(1)/128, mpf('4.7'))
    for k in range(KK+1):
        S_B128[k] += w*vB[k]
S_B128 = [48*v for v in S_B128]
print(f"region B inner 1/128 done, {time.time()-t0:.0f}s")

J = json.load(open("tetra_SB_inner6496.json"))
SB64 = [mpf(s) for s in J["SB_inner64"]]; SB96 = [mpf(s) for s in J["SB_inner96"]]
for k in range(5):
    print(f"k={k}: B128 = {mp.nstr(S_B128[k],58)}")
    print(f"     |B128-B96| = {mp.nstr(abs(S_B128[k]-SB96[k]),3)}   |B96-B64| = {mp.nstr(abs(SB96[k]-SB64[k]),3)}")
J["SB_inner128"] = [mp.nstr(v,70) for v in S_B128]
json.dump(J, open("tetra_SB_inner6496.json","w"), indent=1)
print("saved")
