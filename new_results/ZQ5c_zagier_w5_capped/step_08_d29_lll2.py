import sys, json, time, pickle
from fractions import Fraction as Ff
import mpmath as mp
import flint
sys.path.insert(0,'.')
import w5blib as W
J = json.load(open('out/d29/candidates.json'))
import pickle
A = pickle.load(open('out/d29/atoms5.pkl','rb'))
atom_names=A['atom_names']; keepc=A['keepc']
mp.mp.dps=2050
atom_vals=[mp.mpf(v) for v in A['atom_vals']]
n=len(atom_vals)
print(f"LLL dim {n}", flush=True)
def run(P):
    mp.mp.dps=P+30
    Cc=mp.power(10,P)
    B=[]
    for i,v in enumerate(atom_vals):
        row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*Cc)); B.append(row)
    t1=time.time(); R=flint.fmpz_mat(B).lll(gram='exact'); dt=time.time()-t1
    mp.mp.dps=2050
    rels=[]
    for i in range(n):
        a=[int(R[i,j]) for j in range(n)]
        if not any(a) or a[0]==0: continue
        res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
        if abs(res)<mp.mpf('1e-900'):
            rels.append((max(abs(x) for x in a), a, res))
    rels.sort(key=lambda t:t[0])
    print(f"P={P}: LLL {dt:.0f}s, rels(resid<1e-900): {len(rels)}", flush=True)
    for mx,a,rs in rels[:5]:
        print(f"   maxcoef {mx} T-coef {a[0]} z5-coef {a[1]} resid {mp.nstr(rs,3)} nterms {sum(1 for x in a if x)}", flush=True)
    out=[]
    for mx,a,rs in rels[:10]:
        out.append(dict(maxcoef=mx, a=[str(x) for x in a], resid=mp.nstr(rs,5)))
    json.dump(dict(P=P, dt=dt, nrel=len(rels), rels=out), open(f'out/d29/lll_P{P}.json','w'))
    return dt, [(a) for mx,a,rs in rels]
dt15, r15 = run(1500)
print(f"[time] P=1500: {dt15:.0f}s", flush=True)
dt18, r18 = run(1800)
print(f"[time] P=1800: {dt18:.0f}s", flush=True)
# common relations (same a-vector up to sign)
def norm(a):
    a=tuple(a)
    if any(a):
        for x in a:
            if x: 
                if x<0: return tuple(-y for y in a)
                return a
    return a
s15={norm(a) for a in r15}; s18={norm(a) for a in r18}
common = s15 & s18
print(f"relations at 1500: {len(r15)}, at 1800: {len(r18)}, common: {len(common)}", flush=True)
with open('progress.log','a') as f:
    f.write(f"step_08_d29_lll: P1500 {len(r15)} rels, P1800 {len(r18)} rels, common {len(common)}; dt {dt15:.0f}/{dt18:.0f}s\n")
