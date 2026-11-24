import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)
vals["T"] = vals["L"]
t0=time.time()
R = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=False)
print(f"D=17 target=L: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
for mx,a,res in R["exact_rels"][:2]:
    print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_58_d17_Ltarget.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_58_d17_Ltarget.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-800:])
with open('progress.log','a') as f: f.write("step_58_d17_Ltarget\n")
