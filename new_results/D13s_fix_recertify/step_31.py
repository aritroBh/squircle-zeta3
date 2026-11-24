
verify = r'''
import math, sys
import mpmath as mp
mp.mp.dps = 50
sys.path.insert(0, '.')
from flint import arb, acb, ctx
ctx.dps = 60
from step_38 import iv, pi, rt2, gF, G1bound, h3encl

# Replicate the fixed cell()'s A1m box over (y,sig) box and test containment of true A1
def A1_box(y1, y2, s1, s2):
    y = iv(y1, y2); sig = iv(s1, s2)
    rho = 1 + y*y*sig
    zz = rt2*rho*y/(2+rt2*rho*y)
    zz_lo = arb(zz.lower()); zz_up = arb(zz.upper())
    fA1_lo = arb((1 - zz_lo.atan()/zz_lo).lower())
    fA1_up = arb((1 - zz_up.atan()/zz_up).upper())
    return fA1_lo, fA1_up, arb(y.lower()), arb(y.upper())

msqrt2 = mp.sqrt(2); mpi = mp.pi
def A1_true(yv, sigv):
    rho = 1 + yv*yv*sigv
    zz = msqrt2*rho*yv/(2+msqrt2*rho*yv)
    return (1 - mp.atan(zz)/zz)/yv

worst = mp.mpf('-9'); npts = 0
cells = [(100.,150.),(150.,225.),(225.,340.),(340.,400.),(400.,600.),(1e4,1.5e4),(1e6,1.5e6),(1e9,1.5e9)]
for (n1,n2) in cells:
    y1 = math.sqrt(math.pi/(2*n2)); y2 = math.sqrt(math.pi/(2*n1))
    NS = 4
    for k in range(NS):
        s1, s2 = k/NS, (k+1)/NS
        flo, fup, ylo, yup = A1_box(y1*(1-1e-12), y2*(1+1e-12), s1, s2)
        blo = float(flo)/float(yup); bhi = float(fup)/float(ylo)
        for iy in range(21):
            yv = y1 + (y2-y1)*iy/20
            for isg in range(21):
                sgv = s1 + (s2-s1)*isg/20
                a1 = A1_true(yv, sgv)
                npts += 1
                viol = max(blo - a1, a1 - bhi)
                worst = max(worst, viol)
print(f"A1 containment over {npts} sample points in 8 cells x 4 sigma-subcells: worst violation = {mp.nstr(worst,8)} (<=0 => PROVED box OK)")
'''
import os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
PY='<home>)/.venv/bin/python'
open(os.path.join(FIX,'verify_A1.py'),'w').write(verify)
import subprocess
r = subprocess.run([PY, 'verify_A1.py'], cwd=FIX, capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-1200:] if r.stderr else "")
