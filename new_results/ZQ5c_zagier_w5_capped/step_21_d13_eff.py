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
raw = pickle.load(open('out/d13/vals_cap_2000.pkl','rb'))
mp.mp.dps=2050
vals={k: mp.mpf(v) for k,v in raw.items()}
li5v=[vals[f'li5:{i}'] for i in range(N)]
lg=[vals[f'log:{names[i]}'] for i in range(d)]
def cv(row, comp):
    v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
    for ms_s, cf_s in comp.items():
        ms=eval(ms_s); cf=Ff(cf_s); pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        v -= mp.mpf(cf.numerator)/cf.denominator*pv
    return v
combos=[cv(r,c) for r,c in zip(K,comps)]
keepc=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1900')]
print(f"d13-cap: kernel {len(K)}, EFFECTIVE nonzero combos {len(keepc)}", flush=True)
mags=sorted(abs(v) for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1900'))
print("combo magnitudes (min/median/max):", mp.nstr(mags[0],3), mp.nstr(mags[len(mags)//2],3), mp.nstr(mags[-1],3), flush=True)
