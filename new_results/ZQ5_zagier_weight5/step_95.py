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
N = len(cands)
mp.mp.dps = 1180
li5v = [vals[f"li5:{i}"] for i in range(N)]
combo_vals = [sum(mp.mpf(c)*v for c,v in zip(row, li5v)) for row in K]
nz = [(j, mp.nstr(abs(v),2)) for j,v in enumerate(combo_vals)]
ez = [j for j,v in enumerate(combo_vals) if abs(v) < mp.mpf('1e-1000')]
print("combos:", len(K), " exact-zero combos:", ez)
keep = [j for j in range(len(K)) if j not in ez]
print("kept combos:", len(keep), " -> LLL dim ~", 2+len(keep)+80)
t0=time.time()
R = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=False)
print(f"P=520: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
for mx,a,res in R["exact_rels"][:3]:
    print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_55b_d17_lll520.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_55b_d17_lll520.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_55b_d17_lll520\n")
