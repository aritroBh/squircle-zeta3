
import time
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_04b.py: driver v2 (outer splits at breakpoint-crossings) VALIDATED: octa k0..3 err 4.9e-62/9.1e-61/1.2e-58/8.0e-57 vs exact; cubo k0 4.8e-60 vs K4b cert, k1..4 ~1e-28..1e-31 (inner log-singularity plateau, same as K4b region B). Needs inner refinement for k>=1.\n")
# quick test: cubo with inner 1/64
test = r'''"""step_04c.py -- cubo inner-refinement trend test."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
import importlib
import step_04b_drv as drv
'''
# Instead of editing step_04b, create a reusable driver module
driver = r'''"""k4c_driver.py -- generic chamber driver with adaptive inner h for breakpoint subintervals."""
import mpmath as mp
from mpmath import mpf
from k4_engine import FastModel
from k4_quad import tanh_sinh_nodes

WMIN = mpf('1e-66'); LMIN = mpf('1e-58')
PI = mp.pi

def outer_splits(normals):
    m = len(normals)
    cands = set()
    def add(t):
        if t > 0:
            ph = mp.atan(t)
            if 0 < ph < PI/4: cands.add(ph)
    ns = [tuple(mpf(x) for x in v) for v in normals]
    for j in range(m):
        nx,ny,nz = ns[j]
        if ny+nz != 0: add(-nx/(ny+nz))
        if ny != 0: add(-nx/ny)
        for k in range(j+1,m):
            mx,my,mz = ns[k]
            A = mz*nx - nz*mx; B = nz*my - mz*ny
            if B != 0: add(A/B)
    return sorted(cands)

def chamber_S(normals, K, HO, HI, TO, TI, HI2=None):
    """HI2: finer inner h used for subintervals touching an internal breakpoint."""
    if HI2 is None: HI2 = HI
    M = FastModel(normals, K)
    m = len(normals)
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len':0, 'nbp':0, 'nsing':0}
    splits = [mpf(0)] + outer_splits(normals) + [PI/4]
    for si in range(len(splits)-1):
        phi0, phi1 = splits[si], splits[si+1]
        for (ph, w) in tanh_sinh_nodes(phi0, phi1, HO, TO, WMIN):
            cph, sph = mp.cos(ph), mp.sin(ph)
            if sph <= 0: continue
            thlo = mp.atan(1/sph)
            if thlo >= PI/2 - LMIN:
                stats['zero_len'] += 1
                continue
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
            merged = []
            for b in bps:
                if merged and b - merged[-1] <= LMIN: continue
                merged.append(b)
            bounds = [thlo] + merged + [PI/2]
            stats['nbp'] += len(merged)
            for i in range(len(bounds)-1):
                lo, hi = bounds[i], bounds[i+1]
                if hi - lo <= LMIN:
                    stats['zero_len'] += 1
                    continue
                singular = (i > 0 and i < len(bounds)-2) or (len(merged)>0)
                huse = HI2 if (len(merged)>0) else HI
                if len(merged)>0: stats['nsing'] += 1
                mid = (lo+hi)/2
                sth_, cth_ = mp.sin(mid), mp.cos(mid)
                ux_, uy_, uz_ = sth_*cph, sth_*sph, cth_
                sig = tuple(1 if M.n[j][0]*ux_+M.n[j][1]*uy_+M.n[j][2]*uz_ > 0 else -1 for j in range(m))
                def g(th):
                    sth, cth = mp.sin(th), mp.cos(th)
                    ux, uy, uz = sth*cph, sth*sph, cth
                    Fs = M.F_series(ux, uy, uz, sig)
                    return [Fs[k]*sth for k in range(K+1)]
                for (x, w2) in tanh_sinh_nodes(lo, hi, huse, TI, WMIN):
                    gx = g(x)
                    for k in range(K+1):
                        acc[k] += w*w2*gx[k]
    return [48*v for v in acc], stats
'''
with open('k4c_driver.py','w') as f: f.write(driver)

test = r'''"""step_04c.py -- cubo inner-refinement trend."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
cert = [mpf(s) for s in J["Sdk"]]
for HI2 in [mpf(1)/64, mpf(1)/96]:
    t0=time.time()
    SC, st = chamber_S(NC, 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'), HI2=HI2)
    el = time.time()-t0
    print(f"inner HI2=1/{int(1/HI2)} ({el:.0f}s) stats={st}")
    for k in range(5):
        print(f"   [d^{k}] err = {mp.nstr(abs(SC[k]-cert[k]),3)}", flush=True)
'''
with open('step_04c.py','w') as f: f.write(test)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_04c.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
