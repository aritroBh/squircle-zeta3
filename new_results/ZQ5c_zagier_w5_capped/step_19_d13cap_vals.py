import sys, json, time, pickle
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
