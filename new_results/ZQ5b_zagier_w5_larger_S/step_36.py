import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, os, json, time
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
from fractions import Fraction

os.makedirs('out', exist_ok=True)
CONFIGS = {
    13: dict(split=[3], inert=[2], B=8),
    29: dict(split=[5], inert=[2], B=8),
    21: dict(split=[5], inert=[2], B=6),
    17: dict(split=[2,13], inert=[], B=6),
}
def qd_to_json(x): return [str(x.a), str(x.b)]
for D, cfg in CONFIGS.items():
    F = W.make_field(D, split=cfg['split'], inert=cfg['inert'])
    t0=time.time()
    cands = W.enum_exceptional_fast(F, B=cfg['B'], time_limit=1500)
    d = len(F['letters'])
    data = [dict(x=qd_to_json(c['x']), vx=c['vx'], v1m=c['v1m'], h=c['h'], xf=float(c['x'].emb(40))) for c in cands]
    os.makedirs(f'out/d{D}', exist_ok=True)
    json.dump(dict(D=D, B=cfg['B'], split=cfg['split'], inert=cfg['inert'],
                   names=F['names'], disp=[W.fmt_qd(g) for g in F['letters']],
                   cands=data), open(f'out/d{D}/candidates.json','w'), indent=1)
    import collections
    print(f"D={D} S ram{F['ram']}+inert{cfg['inert']}+split{cfg['split']} d={d} B={cfg['B']}: {len(cands)} candidates ({time.time()-t0:.0f}s), h dist {dict(sorted(collections.Counter(c['h'] for c in cands).items()))}")
'''
open('step_03_enum.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_03_enum.py'], capture_output=True, text=True, cwd='.', timeout=3000)
print(r.stdout); print("ERR:", r.stderr[-800:])
