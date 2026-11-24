import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json
import mpmath as mp
import w5lib as W
import pipe
D=8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
vals = pipe.load_vals(D, 1150)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
mp.mp.dps = 1200
li5v = [vals[f"li5:{i}"] for i in range(len(cands))]
prods = pipe.prod_values(C, vals)
vmap = {"T": vals["T"], "z5": vals["z5"]}
for i,cd in enumerate(cands): vmap[f"Li5({pipe.fmt_qd(cd['x'])})"] = li5v[i]
for kind,msk,nm,pv in prods: vmap[nm] = pv
tot = mp.mpf(0)
for nm, cf in zip(names, c):
    if cf: tot += cf*vmap[nm]
print("max-5 relation residual at 1200 digits:", mp.nstr(tot, 4))
# expected accident scale for dim 30 at P=830:
print("accident scale ~ C^(-29/30)*2^(30/4) =", mp.nstr(mp.mpf(10)**(-830*29/30)*2**7.5, 3))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',src],capture_output=True,text=True,timeout=900)
print(r.stdout); print("STDERR:", r.stderr[-800:])
with open('progress.log','a') as f: f.write("step_17_spurious: tested tiny relation residual\n")
