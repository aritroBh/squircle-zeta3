import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
J = json.load(open('out/d13/candidates_cap.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
names=J['names']; d=len(names); N=len(cands)
Kk = json.load(open('out/d13/kernel5_cap.json'))
K=Kk['K']; comps=Kk['comps']
keepc = [0,1,2,3,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]
P=2100
mp.mp.dps=P+40
F=W.make_field(13, split=[3], inert=[2])
lg=[mp.log(abs(g.emb(P+40))) for g in F['letters']]
z5=mp.zeta(5); pi2=mp.pi**2
T5=W.Lchi_hurwitz(13,5,P)/mp.sqrt(13)
li5v=[W.li5(c['x'].emb(P+40),P+20) for c in cands]
def cv(row, comp):
    v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
    for ms_s, cf_s in comp.items():
        ms=eval(ms_s); cf=Ff(cf_s); pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        v -= mp.mpf(cf.numerator)/cf.denominator*pv
    return v
combos=[cv(K[j],comps[j]) for j in keepc]
# rebuild atom values (d17 types)
import itertools
atoms=[]
for ms in itertools.combinations_with_replacement(range(d),3):
    v=pi2
    for i in ms: v*=lg[i]
    atoms.append(v)
for i in range(d): atoms.append(pi2**2*lg[i])
D5=13**5
atom_vals=[D5*T5, z5]+combos+atoms
# load relations
for fn,tag in (('out/d13/lll_capd17_P1500.json','1500'),('out/d13/lll_capd17_P1800.json','1800')):
    Jr=json.load(open(fn))
    for ri in range(2):
        a=[int(x) for x in Jr['rels'][ri]['a']]
        res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
        mx=max(abs(x) for x in a)
        print(f"{tag} rel{ri}: T-coef {a[0]} maxcoef {mx} FRESH resid @2100: {mp.nstr(res,4)}", flush=True)
with open('progress.log','a') as f:
    f.write("step_23_d13_freshverify: shortest cap-d17 relations re-evaluated @2100d\n")
