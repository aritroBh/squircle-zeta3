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
P = 1300
mp.mp.dps = P+30
print("per-atom |cached - fresh|:")
fresh = {}
fresh["T"] = W.Lchi_hurwitz(D, P)/mp.sqrt(2)
fresh["z5"] = mp.zeta(5)
for i, cd in enumerate(cands):
    fresh[f"Li5({pipe.fmt_qd(cd['x'])})"] = W.li5_series(cd['x'].emb(P+30), P)
prods = pipe.prod_values(C, vals)
vmap = {"T": vals["T"], "z5": vals["z5"]}
for i, cd in enumerate(cands): vmap[f"Li5({pipe.fmt_qd(cd['x'])})"] = vals[f"li5:{i}"]
for kind,msk,nm,pv in prods:
    vmap[nm] = pv
    # fresh product at P
    mp.mp.dps = P+30
    if kind=="log5":
        fv = mp.mpf(1)
        for ii in msk: fv *= mp.log(C['letters'][ii].emb(P+30))
    elif kind=="pi2log3":
        fv = mp.pi**2
        for ii in msk: fv *= mp.log(C['letters'][ii].emb(P+30))
    else:
        fv = mp.pi**4 * mp.log(C['letters'][msk[0]].emb(P+30))
    fresh[nm] = fv
tot_cached = mp.mpf(0); tot_fresh = mp.mpf(0)
for nm, cf in zip(names, c):
    if not cf: continue
    dv = abs(vmap[nm] - fresh[nm])
    if dv > 0: print(f"  {nm:>26s} |cached-fresh| = {mp.nstr(dv,3)}")
    tot_cached += cf*vmap[nm]; tot_fresh += cf*fresh[nm]
print("cached combo:", mp.nstr(tot_cached, 6))
print("fresh  combo:", mp.nstr(tot_fresh, 6))
# double check T alone
print("T cached:", mp.nstr(vmap['T'], 30))
print("T fresh :", mp.nstr(fresh['T'], 30))
'''
open('steps/step_19_compare.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_19_compare.py'],capture_output=True,text=True,timeout=1800)
print(r.stdout); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_19_compare: cached vs fresh per-atom\n")
