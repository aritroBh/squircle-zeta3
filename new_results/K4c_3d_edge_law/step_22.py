
dbg = r'''"""debug cubo chamber k=0 with midpoint rule."""
import os
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 40
from k4_engine import FastModel
PI = mp.pi
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
M = FastModel(NC, 0)
NPH = 160; NTH = 160
tot = mpf(0); tot_nosplit = mpf(0)
negcount = 0
for i in range(NPH):
    ph = PI/4*(i+mpf('0.5'))/NPH
    cph, sph = mp.cos(ph), mp.sin(ph)
    thlo = mp.atan(1/sph)
    # breakpoint for b4
    bps = []
    for j in range(4):
        nx,ny,nz = M.n[j]
        den = nx*cph+ny*sph
        if den == 0: continue
        t = -nz/den
        if t > 0:
            th = mp.atan(t)
            if thlo < th < PI/2: bps.append(th)
    bounds = sorted([thlo]+bps+[PI/2])
    for bi in range(len(bounds)-1):
        lo, hi = bounds[bi], bounds[bi+1]
        mid = (lo+hi)/2
        sthm, cthm = mp.sin(mid), mp.cos(mid)
        uxm,uym,uzm = sthm*cph, sthm*sph, cthm
        sig = [1 if M.n[j][0]*uxm+M.n[j][1]*uym+M.n[j][2]*uzm>0 else -1 for j in range(4)]
        for kk in range(NTH):
            th = lo + (hi-lo)*(kk+mpf('0.5'))/NTH
            sth, cth = mp.sin(th), mp.cos(th)
            ux,uy,uz = sth*cph, sth*sph, cth
            a = [sig[j]*(M.n[j][0]*ux+M.n[j][1]*uy+M.n[j][2]*uz) for j in range(4)]
            if min(a) <= 0: negcount += 1
            F = M.F_series(ux,uy,uz,sig)[0]
            tot += F*sth*(hi-lo)/NTH
    tot_nosplit += 0
tot *= PI/4/NPH
print("midpoint k=0 *48:", mp.nstr(48*tot,30))
print("target (9+3sqrt3)/4:", mp.nstr((9+3*mp.sqrt(3))/4,30))
print("neg a_j cells:", negcount)
# also check sign pattern counts: how many subintervals have b4 negative vs positive
'''
with open('dbg.py','w') as f: f.write(dbg)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','dbg.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else "exit 0")
