import os, json, time
os.chdir("research/08_theory/K4b_cubo_finish")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
import importlib, k4_engine, k4_quad
importlib.reload(k4_engine); importlib.reload(k4_quad)
from k4_engine import *
from k4_quad import *

# ---------- assemble final coefficients ----------
J48 = json.load(open("tetra_Sdk_h48.json"))
JB  = json.load(open("tetra_SB_inner6496.json"))
S_A = [mpf(s) for s in J48["S_A"]]          # err ~ 1e-62
S_B = [mpf(s) for s in JB["SB_inner128"]]   # err ~ 1e-50
S_tot = [S_A[k]+S_B[k] for k in range(5)]
S1exact = (9+3*mp.sqrt(3))/4
print("FINAL assembled coefficients (err: k=0 ~1e-67, k=1,2 ~5e-52, k=3 ~6e-50, k=4 ~5e-49):")
for k in range(5):
    print(f"[d^{k}]S = {mp.nstr(S_tot[k],58)}")
print("k=0 err vs (9+3√3)/4:", mp.nstr(abs(S_tot[0]-S1exact),3))
json.dump({"Sdk":[mp.nstr(v,60) for v in S_tot],
           "err_bound":["1e-67","5e-52","5e-52","1e-50","5e-49"],
           "provenance":"S_A: h48 both shells; S_B: outer 1/32, inner 1/128; dps 90"},
          open("tetra_Sdk_final.json","w"), indent=1)
print("saved tetra_Sdk_final.json")

# ---------- FD convergence probe for F_direct ----------
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

for p in [mpf('1.005'), mpf('0.995')]:
    t0=time.time()
    v32 = S_direct(p, mpf(1)/32, mpf(1)/32)
    t1=time.time()
    v48 = S_direct(p, mpf(1)/32, mpf(1)/48)
    t2=time.time()
    print(f"p={p}: S(1/32)={mp.nstr(v32,50)}  ({t1-t0:.0f}s)")
    print(f"        S(1/48)={mp.nstr(v48,50)}  ({t2-t1:.0f}s)   diff={mp.nstr(abs(v48-v32),3)}")
