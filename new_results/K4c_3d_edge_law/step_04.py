"""step_04.py -- generic O_h-chamber driver; validate on octa & cubo [d^k], k=0..4."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4_engine import FastModel
from k4_quad import tanh_sinh_nodes

WMIN = mpf('1e-66'); LMIN = mpf('1e-58')
PI = mp.pi

def chamber_S(normals, K, HO, HI, TO, TI):
    """Surface series coefficients [d^k]S, k=0..K, using chamber x>=y>=z>=0 (mult 48).
    phi in [0,pi/4], th in [th_lo(phi), pi/2]; generic breakpoints from each normal."""
    M = FastModel(normals, K)
    m = len(normals)
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len':0, 'nbp':0}
    for (ph, w) in tanh_sinh_nodes(mpf(0), PI/4, HO, TO, WMIN):
        cph, sph = mp.cos(ph), mp.sin(ph)
        if sph <= 0: continue
        thlo = mp.atan(1/sph)
        if thlo >= PI/2 - LMIN:
            stats['zero_len'] += 1
            continue
        # breakpoints where <n_j,u>=0 inside (thlo, pi/2)
        bps = []
        for j in range(m):
            nx, ny, nz = M.n[j]
            den = nx*cph + ny*sph
            if den == 0: continue
            t = -nz/den
            if t > 0:
                th = mp.atan(t)
                if th > thlo + LMIN and th < PI/2 - LMIN:
                    bps.append(th)
        bps.sort()
        # merge near-duplicates
        merged = []
        for b in bps:
            if merged and b - merged[-1] <= LMIN:
                continue
            merged.append(b)
        bounds = [thlo] + merged + [PI/2]
        stats['nbp'] += len(merged)
        for i in range(len(bounds)-1):
            lo, hi = bounds[i], bounds[i+1]
            if hi - lo <= LMIN:
                stats['zero_len'] += 1
                continue
            mid = (lo+hi)/2
            # sign pattern at mid
            sth_, cth_ = mp.sin(mid), mp.cos(mid)
            ux_, uy_, uz_ = sth_*cph, sth_*sph, cth_
            sig = []
            for j in range(m):
                nx, ny, nz = M.n[j]
                sig.append(1 if nx*ux_+ny*uy_+nz*uz_ > 0 else -1)
            def g(th):
                sth, cth = mp.sin(th), mp.cos(th)
                ux, uy, uz = sth*cph, sth*sph, cth
                Fs = M.F_series(ux, uy, uz, sig)
                return [Fs[k]*sth for k in range(K+1)]
            for (x, w2) in tanh_sinh_nodes(lo, hi, HI, TI, WMIN):
                gx = g(x)
                for k in range(K+1):
                    acc[k] += w*w2*gx[k]
    return [48*v for v in acc], stats

# ---------------- OCTAHEDRON ----------------
t0 = time.time()
NO = [[1,0,0],[0,1,0],[0,0,1]]
SO, st = chamber_S(NO, 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'))
sq3 = mp.sqrt(3)
exact_O = [4*sq3, 4*sq3, 4*sq3*(1-mp.pi**2/9), 4*sq3*(1+mp.pi**2/9-16*mp.zeta(3)/9), None]
print(f"OCTA h=1/32 ({time.time()-t0:.0f}s) stats={st}")
for k in range(4):
    print(f"  [d^{k}]S = {mp.nstr(SO[k],45)}")
    print(f"     err vs exact = {mp.nstr(abs(SO[k]-exact_O[k]),3)}")

# ---------------- CUBOCTAHEDRON ----------------
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
t0 = time.time()
SC, st = chamber_S(NC, 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'))
print(f"CUBO h=1/32 ({time.time()-t0:.0f}s) stats={st}")
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
cert = [mpf(s) for s in J["Sdk"]]
for k in range(5):
    d = abs(SC[k]-cert[k])
    print(f"  [d^{k}]S = {mp.nstr(SC[k],45)}   |this - K4b cert| = {mp.nstr(d,3)}")
