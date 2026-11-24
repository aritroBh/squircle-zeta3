import os
cwd = os.getcwd()
# Update summary.md with D=29 findings
summary = open('summary.md').read()
summary += """
## D=29 weight-5: NONE FOUND (structural)
- S={29,2,5} (task pool, out/d29/candidates.json, 77 cands, B=8; PROVEN exhausted: B=10,B=12 give same 77):
  kernel dim 7, only 1 EFFECTIVE nonzero combo (6 exact-zero). LLL dim 75 (D17 atom types pi2*log3,pi4*log).
  P=1500 (76s) & P=1800 (95s): all spurious (maxcoef~1e21/1e27), NO common relation. Reach 1500/75=20, 1800/75=24 digits.
- S={29,2,5,7} (enlarged, out/d29/candidates_S7_B5.json, 1088 cands B=5):
  kernel dims h<=2:44, h<=4:78, full:106, but EFFECTIVE nonzero combos = 1 in all cases (K0 ~ -0.507067).
  => zero-symbol Li5 span is only 1-dimensional for D=29. T would have to be c0*z5 + c1*K0 + products.
  h<=2 pool, D17 atom types: dim 94, P=1500(161s)/1800(220s): none, reach 16/19.
  h<=2 pool, FULL atom list (incl z3/T3 products): dim 152, P=1500(757s)/1800(901s): none, reach 10/12.
- CONCLUSION (NUMERICAL): none found for D=29 with S={29,2,5} at 1800 digits (coefficients <~1e19);
  none found for D=29 with S={29,2,5,7} at 1800 digits (coefficients <~1e12).
  Root cause: 1-dimensional effective Li5 span (vs D=13 full pool's 211 effective combos).

## D=13 weight-5: capped pool
- out/d13/candidates_cap.json: 275 cands (grown from 150 to reach kernel>=50), kernel dim 51, rank 224. Values next.
"""
open('summary.md','w').write(summary)
with open('progress.log','a') as f:
    f.write("step_18_done: D=29 NONE (S={29,2,5} and S={29,2,5,7}); 1-dim effective Li5 span; recorded in summary.md\n")
print("summary updated; now launching D=13 capped values")

# Launch D=13 capped values at 2000 digits
vscript = r'''import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
J = json.load(open('out/d13/candidates_cap.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
P=2000
mp.mp.dps=P+40
t0=time.time()
F=W.make_field(13, split=[3], inert=[2])
lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
vals={}
vals['pi2']=mp.pi**2; vals['z5']=mp.zeta(5); vals['z3']=mp.zeta(3)
vals['T']=W.Lchi_hurwitz(13,5,P)/mp.sqrt(13)
vals['L']=W.Lchi_hurwitz(13,5,P)
vals['T3']=W.Lchi_hurwitz(13,3,P)/mp.sqrt(13)
for i,l in enumerate(lg): vals[f'log:{names[i]}']=l
print(f"base vals done ({time.time()-t0:.0f}s); L={mp.nstr(vals['L'],20)}", flush=True)
li5v=[None]*N
for i,c in enumerate(cands):
    x=c['x'].emb(P+40)
    v=W.li5(x,P+20)
    li5v[i]=v; vals[f'li5:{i}']=v
    if (i+1)%50==0: print(f"  li5 {i+1}/{N} ({time.time()-t0:.0f}s)", flush=True)
print(f"values @2000d: {N} Li5 done ({time.time()-t0:.0f}s)", flush=True)
pickle.dump({k: mp.nstr(v,P+35) for k,v in vals.items()}, open('out/d13/vals_cap_2000.pkl','wb'))
print("saved out/d13/vals_cap_2000.pkl", flush=True)
# cross-check
mp.mp.dps=P+190
worst=mp.mpf(0); nw=0
for i,c in enumerate(cands):
    x=c['x'].emb(P+190)
    ax=abs(x)
    if ax<=mp.mpf('0.9'): v2=W.li5_series(x, mp.mpf(10)**(-(P+150)))
    else: v2=W.li5_log(x, mp.mpf(10)**(-(P+150)))
    mp.mp.dps=P+40
    dd=abs(v2-li5v[i]); mp.mp.dps=P+190
    if dd>worst: worst=dd; nw=i
print(f"Li5 cross-check: worst {mp.nstr(worst,3)} at {nw}", flush=True)
xmap={}
for i,c in enumerate(cands): xmap[c['x']]=i
ndup=0; wd=mp.mpf(0)
mp.mp.dps=P+40
for i,c in enumerate(cands):
    x2=c['x']*c['x']; xm=-c['x']
    if x2 in xmap and xm in xmap:
        ndup+=1
        r=li5v[xmap[x2]]-16*(li5v[i]+li5v[xmap[xm]])
        if abs(r)>wd: wd=abs(r)
print(f"duplication pairs: {ndup}, worst {mp.nstr(wd,3)}", flush=True)
with open('progress.log','a') as f:
    f.write(f"step_19_d13cap_vals: values @2000d (275 Li5); 2nd-way worst {mp.nstr(worst,2)}; dup {ndup} worst {mp.nstr(wd,2)}\n")
print("D13 VALS DONE", flush=True)
'''
open('step_19_d13cap_vals.py','w').write(vscript)
import subprocess, time
logf=open('step_19_bg.log','w')
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','step_19_d13cap_vals.py'],
                      stdout=logf, stderr=subprocess.STDOUT, cwd=cwd)
open('step_19_pid.txt','w').write(str(proc.pid))
print("launched D=13 values PID", proc.pid, time.strftime('%H:%M:%S'))
