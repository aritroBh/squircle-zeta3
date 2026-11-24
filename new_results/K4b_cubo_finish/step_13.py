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

def make_driver(HO, HI, TO, TI):
    stats = {'zero_len':0, 'badnode':0}
    def region_vec(ph, region):
        thB = th_B(ph); thA = th_A(ph)
        cph, sph = mp.cos(ph), mp.sin(ph)
        if region == 'B':
            lo, hi, sig = mpf(0), thB, SIG_B
        else:
            lo, hi, sig = thB, thA, SIG_A
        if hi - lo <= LMIN:
            stats['zero_len'] += 1
            return [mpf(0)]*(KK+1)
        def g(th):
            sth, cth = mp.sin(th), mp.cos(th)
            ux, uy, uz = sth*cph, sth*sph, cth
            for j in range(4):
                aj = sig[j]*(MT.n[j][0]*ux + MT.n[j][1]*uy + MT.n[j][2]*uz)
                if aj <= 0:
                    stats['badnode'] += 1
                    return None
            Fs = MT.F_series(ux, uy, uz, sig)
            return [Fs[k]*sth for k in range(KK+1)]
        acc = [mpf(0)]*(KK+1)
        for (x,w) in tanh_sinh_nodes(lo, hi, HI, TI, WMIN):
            gx = g(x)
            if gx is None: continue
            for k in range(KK+1): acc[k] += w*gx[k]
        return acc
    def run():
        S_A = [mpf(0)]*(KK+1); S_B = [mpf(0)]*(KK+1)
        for (x,w) in tanh_sinh_nodes(mp.pi/4, mp.pi/2, HO, TO, WMIN):
            vB = region_vec(x, 'B'); vA = region_vec(x, 'A')
            for k in range(KK+1):
                S_B[k] += w*vB[k]; S_A[k] += w*vA[k]
        S_A = [48*v for v in S_A]; S_B = [48*v for v in S_B]
        return S_A, S_B, [S_A[k]+S_B[k] for k in range(KK+1)]
    return run, stats

run, stats = make_driver(mpf(1)/48, mpf(1)/48, mpf('4.7'), mpf('4.7'))
t0=time.time()
S_A48, S_B48, S_tot48 = run()
el = time.time()-t0
print(f"h=1/48 run, {el:.0f}s, stats={stats}")
S1exact = (9+3*mp.sqrt(3))/4
print(f"k=0 err vs (9+3√3)/4: {mp.nstr(abs(S_tot48[0]-S1exact),3)}")
# load h=1/32 and compare
prev = json.load(open("tetra_Sdk_h32.json"))
S_tot32 = [mpf(s) for s in prev["S_tot"]]
for k in range(5):
    d = abs(S_tot48[k]-S_tot32[k])
    print(f"[d^{k}] S_tot(h48) = {mp.nstr(S_tot48[k],65)}")
    print(f"      |h48-h32| = {mp.nstr(d,3)}")
json.dump({"h":"1/48","T":"4.7","dps":90,
           "S_A":[mp.nstr(v,70) for v in S_A48],"S_B":[mp.nstr(v,70) for v in S_B48],
           "S_tot":[mp.nstr(v,70) for v in S_tot48],"stats":stats,
           "k0_err": mp.nstr(abs(S_tot48[0]-S1exact),3)},
          open("tetra_Sdk_h48.json","w"), indent=1)
print("saved tetra_Sdk_h48.json")
