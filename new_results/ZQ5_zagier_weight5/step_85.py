import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, os
import mpmath as mp
import w5lib as W
import pipe
for D in (21, 13, 29):
    C = W.letters_info(D)
    J, cands = pipe.load_cands(D)
    kf = pipe.cand_dir(D)+"/kernel_b10.json" if D==21 else pipe.cand_dir(D)+"/kernel.json"
    K = [list(map(int, r)) for r in json.load(open(kf))["K"]]
    if D==21:
        raw = __import__('pickle').load(open(pipe.cand_dir(D)+"/vals_b10.pkl","rb"))
        old = mp.mp.dps; mp.mp.dps = 1250
        vals = {k: mp.mpf(v) for k, v in raw.items()}
        mp.mp.dps = old
    else:
        vals = pipe.load_vals(D, 1150)
    vals["T"] = vals["L"]     # target L itself, no sqrt scaling
    for P in (830,):
        R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
        print(f"D={D} target=L: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])}")
        for mx,a,res in R["exact_rels"][:2]:
            print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_52_Ltarget.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_52_Ltarget.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_52_Ltarget: L-target variants for D=21/13/29\n")
