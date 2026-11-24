import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
src=r'''
"""Weight-12 plane test: for points on Pi1 (s,s,1,1) and Pi2 (s,1,s,1), s=1..15,
compute c12 value coords in the W=12 basis, and the two quotient obstructions:
  obs1 = z39 + (14/9)*z57   (coefficient of z39 in MZV_12/span(products), using z57==(14/9)z39 mod products)
  obs2 = w*
Also membership rank test. Saves to w12_planes.json."""
import json, time
from fractions import Fraction as F
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES:
        E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
out=[]
for plane in ['Pi1','Pi2']:
    for s in range(1,16):
        a=(s,s,1,1) if plane=='Pi1' else (s,1,s,1)
        t=time.time()
        r=E.reduced_coeffs(a,12)[12]
        obs1=r[9]+F(14,9)*r[10]
        out.append({'plane':plane,'s':s,'z39':str(r[9]),'z57':str(r[10]),'w*':str(r[11]),
                    'z2^2z35':str(r[7]),'z2z37':str(r[8]),'obs1':str(obs1),'sec':round(time.time()-t,1)})
        json.dump(out, open('w12_planes.json','w'))
print('done', len(out))
'''
open('w12_planes.py','w').write(src)
import subprocess, sys
proc=subprocess.Popen(['nice','-n','10','nohup',sys.executable,'w12_planes.py','>','w12_planes.out','2>&1'],
                      cwd=WORK, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print('launched', proc.pid)
open('progress.log','a').write('step_59-65: W=12 table (nullity 12, duality, relation z57=(128195/12438)z12-(14/9)z3z9-(25/3)z5z7+(14/9)z39 numeric 6e-61); plane probes: z39,z57 nonzero on Pi1/Pi2/diag BUT ratio -14/9 exact => cancels via relation => still products; w*,z2^2z35,z2z37 zero on loci; 15-pt verification launched\n')
