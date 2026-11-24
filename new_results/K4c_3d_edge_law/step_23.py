
step04b = r'''"""step_04b.py -- generic chamber driver v2 with outer split points; validate octa & cubo."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4_engine import FastModel
from k4_quad import tanh_sinh_nodes

WMIN = mpf('1e-66'); LMIN = mpf('1e-58')
PI = mp.pi

def outer_splits(normals, M=None):
    """phi in (0,pi/4) where a theta-breakpoint enters/leaves [th_lo, pi/2] or two meet."""
    m = len(normals)
    cands = set()
    def add(t):
        if t > 0:
            ph = mp.atan(t)
            if 0 < ph < PI/4: cands.add(ph)
    ns = [tuple(mpf(x) for x in v) for v in normals]
    for j in range(m):
        nx,ny,nz = ns[j]
        # theta_C = th_lo: tan(phi) = -nx/(ny+nz)
        if ny+nz != 0: add(-nx/(ny+nz))
        # theta_C = pi/2: den = nx cos + ny sin = 0: tan(phi) = -nx/ny
        if ny != 0: add(-nx/ny)
        # th_lo itself degenerate when 1/sin -> ... skip
        for k in range(j+1,m):
            mx,my,mz = ns[k]
            # theta_Cj = theta_Ck: (mz*nx - nz*mx) cos = (nz*my - mz*ny) sin
            A = mz*nx - nz*mx; B = nz*my - mz*ny
            if B != 0: add(A/B)
    return sorted(cands)

def chamber_S(normals, K, HO, HI, TO, TI):
    M = FastModel(normals, K)
    m = len(normals)
    acc = [mpf(0)]*(K+1)
    stats = {'zero_len':0, 'nbp':0}
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
                mid = (lo+hi)/2
                sth_, cth_ = mp.sin(mid), mp.cos(mid)
                ux_, uy_, uz_ = sth_*cph, sth_*sph, cth_
                sig = tuple(1 if M.n[j][0]*ux_+M.n[j][1]*uy_+M.n[j][2]*uz_ > 0 else -1 for j in range(m))
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

# OCTA
t0=time.time()
SO, st = chamber_S([[1,0,0],[0,1,0],[0,0,1]], 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'))
sq3 = mp.sqrt(3)
exact_O = [4*sq3, 4*sq3, 4*sq3*(1-mp.pi**2/9), 4*sq3*(1+mp.pi**2/9-16*mp.zeta(3)/9)]
print(f"OCTA ({time.time()-t0:.0f}s) stats={st}")
for k in range(4):
    print(f"  [d^{k}] err = {mp.nstr(abs(SO[k]-exact_O[k]),3)}")

# CUBO
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
print("outer split points for cubo:", [mp.nstr(x,8) for x in outer_splits(NC)], " atan(1/2)=", mp.nstr(mp.atan(mpf('0.5')),8))
t0=time.time()
SC, st = chamber_S(NC, 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'))
print(f"CUBO ({time.time()-t0:.0f}s) stats={st}")
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
cert = [mpf(s) for s in J["Sdk"]]
for k in range(5):
    d = abs(SC[k]-cert[k])
    print(f"  [d^{k}] |this - K4b cert| = {mp.nstr(d,3)}")
'''
with open('step_04b.py','w') as f: f.write(step04b)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_04b.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "exit 0")
