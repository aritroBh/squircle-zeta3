import sys, subprocess
p = open('pipe.py').read()
p = p.replace("    logs = [mp.log(g.emb(P+30)) for g in C['letters']]",
              "    logs = [mp.log(abs(g.emb(P+30))) for g in C['letters']]")
open('pipe.py','w').write(p)
src = r'''
import sys; sys.path.insert(0,'.')
import pickle, time
import mpmath as mp
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
# load old cache carefully: keep only real-valued entries (li5:*, L, T, z5, pi2)
raw = pickle.load(open(pipe.cand_dir(D)+"/vals_1150.pkl","rb"))
mp.mp.dps = 1250
out = {}
for k, v in raw.items():
    if k.startswith("li5:") or k in ("L","T","z5","pi2"):
        out[k] = v   # strings, real
# recompute logs with abs
mp.mp.dps = 1180
for n, g in zip(C['names'], C['letters']):
    out["log:"+n] = mp.nstr(mp.log(abs(g.emb(1180))), 1175)
pickle.dump(out, open(pipe.cand_dir(D)+"/vals_1150.pkl","wb"))
print("cache repaired:", len(out), "entries")
vals = pipe.load_vals(D, 1150)
mp.mp.dps = 60
print("logs:", [mp.nstr(vals["log:"+n], 12) for n in C['names']])
print("T =", mp.nstr(vals["T"], 15))
# quick LLL at 520
import json
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
t0=time.time()
R = pipe.phase_lll(D, C, cands_g := pipe.load_cands(D)[1], K, vals, P=520, use_atoms=False)
print(f"P=520: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
for mx,a,res in R["exact_rels"][:3]:
    print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_56_d17_fix.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_56_d17_fix.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_56_d17_fix: abs log; cache repaired; LLL 520\n")
