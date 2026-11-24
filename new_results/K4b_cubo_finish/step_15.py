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
    lo, hi, sig = mpf(0), thB, SIG_B
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
    for (x,w) in tanh_sinh_nodes(lo, hi, HI, TI, WMIN):
        gx = g(x)
        if gx is None: continue
        for k in range(KK+1): acc[k] += w*gx[k]
    return acc

def run_regionB(HO, TO, HI, TI):
    S_B = [mpf(0)]*(KK+1)
    for (x,w) in tanh_sinh_nodes(mp.pi/4, mp.pi/2, HO, TO, WMIN):
        vB = regionB_vec(x, HI, TI)
        for k in range(KK+1):
            S_B[k] += w*vB[k]
    return [48*v for v in S_B]

res = {}
for HI in [mpf(1)/64, mpf(1)/96]:
    t0=time.time()
    SB = run_regionB(mpf(1)/32, mpf('4.7'), HI, mpf('4.7'))
    res[str(HI)] = SB
    print(f"region B, outer 1/32, inner {HI} ({time.time()-t0:.0f}s):")
    for k in range(5):
        print(f"   k={k}: {mp.nstr(SB[k],60)}")

# compare with previous region B (h32: outer/inner 1/32; h48: 1/48)
J32 = json.load(open("tetra_Sdk_h32.json")); J48 = json.load(open("tetra_Sdk_h48.json"))
SB32 = [mpf(s) for s in J32["S_B"]]; SB48 = [mpf(s) for s in J48["S_B"]]
SB64 = res[str(mpf(1)/64)]; SB96 = res[str(mpf(1)/96)]
print("\ninner-convergence (outer=1/32 for 32,64,96):")
for k in range(5):
    print(f"k={k}: |B64-B32|={mp.nstr(abs(SB64[k]-SB32[k]),3)}  |B96-B64|={mp.nstr(abs(SB96[k]-SB64[k]),3)}  |B48(outer48)-B32|={mp.nstr(abs(SB48[k]-SB32[k]),3)}")
json.dump({"SB_inner64":[mp.nstr(v,70) for v in SB64],
           "SB_inner96":[mp.nstr(v,70) for v in SB96],
           "outer":"1/32","T":"4.7"}, open("tetra_SB_inner6496.json","w"), indent=1)
print("saved tetra_SB_inner6496.json")
