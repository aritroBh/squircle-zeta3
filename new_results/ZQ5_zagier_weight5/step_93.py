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
t0=time.time()
R = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=False)
print(f"P=520: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
for mx,a,res in R["exact_rels"][:3]:
    print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_55_d17_lll520.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_55_d17_lll520.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
with open('progress.log','a') as f: f.write("step_55_d17_lll520\n")
