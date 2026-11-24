import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('progress.log','a') as f:
    f.write('step_06_w3_d5: T3=L(chi5,3)/sqrt5 in Li3 span (T-coef 125, resid 7.4e-525); L-target none => T3 = L/sqrt D\n')
    f.write('step_07_w3_d13: T3=L(chi13,3)/sqrt13=0.253084162506679790977235239912 in Li3 span (T-coef 2197=13^3, resid -1.3e-513); L-target none. Control done\n')
# step_08: enumerate D=13 weight-5 pool (B=8) and save
step = r'''
import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
import os
os.makedirs('out/d13', exist_ok=True)
F = W.make_field(13, split=[3], inert=[2])
t0=time.time()
cands = W.enum_exceptional_fast(F, B=8)
data=[dict(x=[str(c['x'].a),str(c['x'].b)], vx=c['vx'], v1m=c['v1m'], h=c['h'], xf=float(c['x'].emb(40))) for c in cands]
json.dump(dict(D=13, B=8, split=[3], inert=[2], names=F['names'], disp=[str(W.fmt_qd(g)) for g in F['letters']],
               S="S={13,2,3} (13 ram, 2 inert, 3 split: (5+-sqrt13)/2)", cands=data), open('out/d13/candidates.json','w'), indent=1)
print(f"D=13: saved {len(data)} exceptional S-units (B=8) ({time.time()-t0:.0f}s)")
import collections
print("h dist:", dict(sorted(collections.Counter(c['h'] for c in cands).items())))
'''
open('step_08_d13_enum.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_08_d13_enum.py'], capture_output=True, text=True, cwd='.', timeout=1700)
print(r.stdout); print("ERR:", r.stderr[-800:])
