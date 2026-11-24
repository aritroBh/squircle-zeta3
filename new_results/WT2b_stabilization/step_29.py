base='<project>/research/08_theory/WT2_weight2_all_k'
import os
os.chdir(base)
log=open('progress.log','a'); log.write('\n--- RESUMED RUN: recovering state; record2.pkl has W,H,Gm for k<=12; goal: prove Stabilization Lemma for all k ---\n'); log.close()
step = r'''import sys, importlib, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *

record = pickle.load(open('record2.pkl','rb'))

def real_form(Gm):
    return sp.factor(sp.together(pf_to_sympy(Gm)))

# For each canonical key, sequence of Gm over k, and iterated differences (in k)
def pf_is_zero(P):
    return (not P[0]) and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

keys = sorted(set(key for (k,key) in record), key=lambda x:(sum(x), -x[0]))
for key in keys:
    d = sum(key)
    if d < 1 or d > 6: continue
    ks = [k for k in range(d, 13) if (k,key) in record and 'Gm' in record[(k,key)]]
    if len(ks) < 4: continue
    gms = [record[(k,key)]['Gm'] for k in ks]
    # iterated k-differences: find smallest m such that delta^m is eventually zero
    seq = list(gms)
    stab = None
    for m in range(0, 6):
        # check if seq is eventually constant
        tail_diffs = [pf_sub(seq[t+1], seq[t]) for t in range(len(seq)-1)]
        # find last nonzero diff index
        lastnz = -1
        for t, df in enumerate(tail_diffs):
            if not pf_is_zero(df): lastnz = t
        # seq constant from index lastnz+1 on, i.e. for k >= ks[lastnz+1]
        if lastnz <= 1:
            stab = (m, ks[lastnz+1] if lastnz>=0 else ks[0])
            break
        seq = tail_diffs
    print(f"key {key}: (1-eps)^{m} structure with m={stab[0]}, stable for k>={stab[1]}" if stab else f"key {key}: no stab found <=5")
'''
open('step_61.py','w').write(step)
import subprocess
print(subprocess.run(['<project>/.venv/bin/python','step_61.py'],capture_output=True,text=True,cwd=base).stdout)
