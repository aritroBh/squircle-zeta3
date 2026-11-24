import os
cwd = os.getcwd()
escript = r'''import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F=W.make_field(29, split=[5], inert=[2])
for B in (10,12):
    t0=time.time()
    c = W.enum_exceptional_fast(F, B=B, time_limit=1200)
    from collections import Counter
    hc = Counter(x['h'] for x in c)
    print(f"B={B}: {len(c)} exceptional, heights {dict(sorted(hc.items()))} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(D=29, B=B, split=[5], inert=[2], names=F['names'],
                   disp=[W.fmt_qd(g) for g in F['letters']], S=[29,2,5],
                   cands=[dict(x=[str(x['x'].a), str(x['x'].b)], vx=x['vx'], v1m=x['v1m'], h=x['h'],
                               xf=float(x['x'].emb(40))) for x in c]),
              open(f'out/d29/candidates_B{B}.json','w'))
    print(f"saved out/d29/candidates_B{B}.json", flush=True)
'''
open('step_09_d29_enum.py','w').write(escript)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_09_d29_enum.py'],
                   capture_output=True, text=True, cwd=cwd, timeout=2700)
print("elapsed", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
