import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
print("start", flush=True)
import json, time
import mpmath as mp
import w5lib as W
import pipe
print("imports ok", flush=True)
D=17
C = W.letters_info(D)
print("letters ok", flush=True)
J, cands = pipe.load_cands(D)
print("cands ok", len(cands), flush=True)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
print("K ok", len(K), flush=True)
vals = pipe.load_vals(D, 1150)
print("vals ok", len(vals), flush=True)
mp.mp.dps = 60
print("T =", mp.nstr(vals["T"], 10), flush=True)
print("li5:0 =", mp.nstr(vals["li5:0"], 10), flush=True)
mp.mp.dps = 1180
li5v = [vals[f"li5:{i}"] for i in range(len(cands))]
print("li5v ok", flush=True)
cv = [sum(mp.mpf(c)*v for c,v in zip(row, li5v)) for row in K]
print("combos ok", len(cv), flush=True)
prods = pipe.prod_values(C, vals)
print("prods ok", len(prods), flush=True)
'''
open('steps/step_55c_d17_trace.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_55c_d17_trace.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
