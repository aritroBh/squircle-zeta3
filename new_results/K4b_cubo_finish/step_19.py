import os, json, time
os.chdir("research/08_theory/K4b_cubo_finish")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 65
import importlib, k4_engine, k4_quad
importlib.reload(k4_engine); importlib.reload(k4_quad)
from k4_engine import *
from k4_quad import *

WMIN = mpf('1e-68'); LMIN = mpf('1e-60')
MT = fast_tetra(K=4)

def S_direct(p, HO, HI, TO=mpf('4.7'), TI=mpf('4.7')):
    d = p - 1
    def region_int(ph, region):
        thB = th_B(ph); thA = th_A(ph)
        cph, sph = mp.cos(ph), mp.sin(ph)
        if region == 'B':
            lo, hi, sig = mpf(0), thB, SIG_B
        else:
            lo, hi, sig = thB, thA, SIG_A
        if hi - lo <= LMIN:
            return mpf(0)
        def g(th):
            sth, cth = mp.sin(th), mp.cos(th)
            ux, uy, uz = sth*cph, sth*sph, cth
            for j in range(4):
                aj = sig[j]*(MT.n[j][0]*ux + MT.n[j][1]*uy + MT.n[j][2]*uz)
                if aj <= 0:
                    return None
            return MT.F_direct(ux, uy, uz, sig, d)*sth
        acc = mpf(0)
        for (x,w) in tanh_sinh_nodes(lo, hi, HI, TI, WMIN):
            gx = g(x)
            if gx is None: continue
            acc += w*gx
        return acc
    s = mpf(0)
    for (x,w) in tanh_sinh_nodes(mp.pi/4, mp.pi/2, HO, TO, WMIN):
        s += w*(region_int(x,'B') + region_int(x,'A'))
    return 48*s

EPS = mpf('0.005')
ORDER = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6]
out = {}
if os.path.exists("fd_Svals.json"):
    out = json.load(open("fd_Svals.json"))
for j in ORDER:
    key = str(j)
    if key in out: continue
    t0 = time.time()
    v = S_direct(1 + j*EPS, mpf(1)/32, mpf(1)/48)
    out[key] = mp.nstr(v, 55)
    json.dump(out, open("fd_Svals.json","w"), indent=1)
    pj = 1 + j*EPS
    print(f"j={j:+d}  p={mp.nstr(pj,8)}  S={mp.nstr(v,45)}  ({time.time()-t0:.0f}s)", flush=True)
print("ALL DONE")
