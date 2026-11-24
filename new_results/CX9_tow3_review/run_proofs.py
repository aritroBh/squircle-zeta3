"""Rerun read-only submitted certificates; all output stays in this folder."""
from pathlib import Path
import os, subprocess, sys, time
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'TOW3_tower_depth6'
cases=[('prove_V56.py',d,s) for d in (5,6) for s in range(d//2+1)]
cases += [('prove_E3b56.py',d) for d in (5,6)]
if len(sys.argv)>1:
    cases=[('prove_adj56.py',int(sys.argv[1]),int(sys.argv[2]))]
for case in cases:
    start=time.monotonic()
    print('RUN',*case,flush=True)
    r=subprocess.run([sys.executable,'-B',str(SRC/case[0]),*map(str,case[1:])],cwd=ROOT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=900)
    assert r.returncode==0,case
    print('PASS',*case,f'{time.monotonic()-start:.2f}s',flush=True)
