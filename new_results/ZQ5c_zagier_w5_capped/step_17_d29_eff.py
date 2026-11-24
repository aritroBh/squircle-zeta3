import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
JF = json.load(open('out/d29/candidates_S7_B5.json'))
full=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in JF['cands']]
names=JF['names']; d=len(names)
raw = pickle.load(open('out/d29/vals_S7_2000.pkl','rb'))
mp.mp.dps=2050
vals={k: mp.mpf(v) for k,v in raw.items()}
lg=[vals[f'log:{names[i]}'] for i in range(d)]
def eff_count(tag):
    Kk = json.load(open(f'out/d29/kernel5_{tag}.json'))
    K=Kk['K']; comps=Kk['comps']; hmax=Kk['hmax']
    idx=[i for i,c in enumerate(full) if c['h']<=hmax]
    li5v=[vals[f'li5:{i}'] for i in idx]
    def cv(row, comp):
        v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
        for ms_s, cf_s in comp.items():
            ms=eval(ms_s); cf=Ff(cf_s)
            pv=mp.mpf(1)
            for ii in ms: pv*=lg[ii]
            v -= mp.mpf(cf.numerator)/cf.denominator*pv
        return v
    combos=[cv(r,c) for r,c in zip(K,comps)]
    nz=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1900')]
    print(f"{tag} (h<={hmax}): kernel {len(K)}, EFFECTIVE nonzero combos {len(nz)}", flush=True)
    # magnitude of nonzero combos
    for j in nz[:5]:
        print(f"   K{j} ~ {mp.nstr(combos[j],6)}", flush=True)
    return len(nz)
eff_count('S7h2')
eff_count('S7h4')
