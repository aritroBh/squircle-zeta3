import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=5
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
for P in (520, 830):
    t0=time.time()
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"P={P}: rels={len(R['rels'])} exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
    if R["exact_rels"]:
        mx,a,res = R["exact_rels"][0]
        print("   best exact: a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_28_d5_lll.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_28_d5_lll.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_28_d5_lll: D=5 numerics + kernel LLL\n")
