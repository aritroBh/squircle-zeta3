"""step_01.py -- K4c: copy K4b engine; octahedron & cuboctahedron k=0 checks."""
import os, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from k4_engine import *
from k4_quad import *

WMIN = mpf('1e-60')
t0 = time.time()

# ---------- Octahedron: 3 coordinate normals, integrate over first octant x8 ----------
MO = fast_octa(K=0)
def octa_k0(h, T):
    s = mpf(0)
    for (ph, w1) in tanh_sinh_nodes(mpf(0), mp.pi/2, h, T, WMIN):
        cph, sph = mp.cos(ph), mp.sin(ph)
        def g(th):
            sth, cth = mp.sin(th), mp.cos(th)
            ux, uy, uz = sth*cph, sth*sph, cth
            Fs = MO.F_series(ux, uy, uz, [1,1,1])
            return Fs[0]*sth
        s += w1*ts_quad(g, mpf(0), mp.pi/2, h, T)
    return 8*s
S_oct = octa_k0(mpf(1)/32, mpf('4.5'))
print(f"OCTA k=0: {mp.nstr(S_oct,40)}")
print(f"   err vs 4*sqrt(3): {mp.nstr(abs(S_oct-4*mp.sqrt(3)),3)}")

# ---------- Cuboctahedron: K4b driver (regions A,B, mult 48), k=0 only ----------
MT = fast_tetra(K=0)
SIG_A = [+1,+1,+1,+1]; SIG_B = [+1,-1,+1,+1]
def th_A(ph): return mp.atan(1/mp.sin(ph))
def th_B(ph): return mp.atan(1/(mp.cos(ph)+mp.sin(ph)))
LMIN = mpf('1e-55')
def cubo_k0(HO, HI, T):
    SA = mpf(0); SB = mpf(0)
    for (ph, w1) in tanh_sinh_nodes(mp.pi/4, mp.pi/2, HO, T, WMIN):
        cph, sph = mp.cos(ph), mp.sin(ph)
        thB, thA = th_B(ph), th_A(ph)
        for (lo, hi, sig) in [(mpf(0), thB, SIG_B), (thB, thA, SIG_A)]:
            if hi-lo <= LMIN: continue
            def g(th):
                sth, cth = mp.sin(th), mp.cos(th)
                ux, uy, uz = sth*cph, sth*sph, cth
                ok = all(sig[j]*(MT.n[j][0]*ux+MT.n[j][1]*uy+MT.n[j][2]*uz) > 0 for j in range(4))
                if not ok: return mpf(0)
                return MT.F_series(ux, uy, uz, sig)[0]*sth
            v = ts_quad(g, lo, hi, HI, T)
            if sig is SIG_A: SA += w1*v
            else: SB += w1*v
    return 48*SA, 48*SB
SA, SB = cubo_k0(mpf(1)/32, mpf(1)/32, mpf('4.5'))
S_cub = SA+SB
S1 = (9+3*mp.sqrt(3))/4
print(f"CUBO k=0: A={mp.nstr(SA,25)}  (vs 3sqrt3/4 err {mp.nstr(abs(SA-3*mp.sqrt(3)/4),3)})")
print(f"          B={mp.nstr(SB,25)}  (vs 9/4 err {mp.nstr(abs(SB-mpf(9)/4),3)})")
print(f"        tot={mp.nstr(S_cub,40)}")
print(f"   err vs (9+3sqrt3)/4: {mp.nstr(abs(S_cub-S1),3)}")
print(f"elapsed {time.time()-t0:.0f}s")
