import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
J, cands = w5pipe.load_cands('out/d29', 29)
names=J['names']; d=len(names); N=len(cands)
P=2000
mp.mp.dps=P+40
t0=time.time()
F=W.make_field(29, split=[5], inert=[2])
lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
vals={}
vals['pi2']=mp.pi**2; vals['z5']=mp.zeta(5); vals['z3']=mp.zeta(3)
vals['T']=W.Lchi_hurwitz(29,5,P)/mp.sqrt(29)
vals['L']=W.Lchi_hurwitz(29,5,P)
vals['T3']=W.Lchi_hurwitz(29,3,P)/mp.sqrt(29)
for i,l in enumerate(lg): vals[f'log:{names[i]}']=l
li5v=[]
for i,c in enumerate(cands):
    x=c['x'].emb(P+40)
    v=W.li5(x,P+20)
    li5v.append(v); vals[f'li5:{i}']=v
print(f"values @2000d: {N} Li5 + L/T/T3 ({time.time()-t0:.0f}s)", flush=True)
print("L(chi29,5) =", mp.nstr(vals['L'],30), flush=True)
print("T = L/sqrt29 =", mp.nstr(vals['T'],30), flush=True)
print("L(chi29,3)/sqrt29 =", mp.nstr(vals['T3'],30), flush=True)
# cross-check every Li5 a second way
mp.mp.dps=P+150+40
worst=mp.mpf(0); nw=0
for i,c in enumerate(cands):
    x=c['x'].emb(P+190)
    ax=abs(x)
    if ax<=mp.mpf('0.9'):
        v2=W.li5_series(x, mp.mpf(10)**(-(P+150)))
    else:
        v2=W.li5_log(x, mp.mpf(10)**(-(P+150)))
    mp.mp.dps=P+40
    dd=abs(v2-li5v[i])
    mp.mp.dps=P+190
    if dd>worst: worst=dd; nw=i
print(f"Li5 cross-check (2nd way): worst diff {mp.nstr(worst,3)} at index {nw}", flush=True)
# duplication relations inside pool
xmap={}
for i,c in enumerate(cands): xmap[c['x']]=i
ndup=0; worst_dup=mp.mpf(0)
mp.mp.dps=P+40
for i,c in enumerate(cands):
    x2=c['x']*c['x']; xm=-c['x']
    if x2 in xmap and xm in xmap:
        ndup+=1
        r=li5v[xmap[x2]]-16*(li5v[i]+li5v[xmap[xm]])
        if abs(r)>worst_dup: worst_dup=abs(r)
print(f"duplication pairs in pool: {ndup}, worst residual {mp.nstr(worst_dup,3)}", flush=True)
pickle.dump({k: mp.nstr(v,P+35) for k,v in vals.items()}, open('out/d29/vals_2000.pkl','wb'))
print("saved out/d29/vals_2000.pkl", flush=True)
with open('progress.log','a') as f:
    f.write(f"step_06_d29_vals: values @2000d; L={mp.nstr(vals['L'],20)}; 2nd-way worst {mp.nstr(worst,2)}; dup pairs {ndup} worst {mp.nstr(worst_dup,2)}\n")
