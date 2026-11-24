
import sys, json, time, pickle, itertools
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
J = json.load(open('out/d13/candidates.json'))
K5 = json.load(open('out/d13/kernel5.json'))
K=K5['K']; comps=K5['comps']; names=K5['names']; d=K5['d']
N=len(J['cands'])
raw = pickle.load(open('out/d13/vals_1400.pkl','rb'))
mp.mp.dps=1450
vals = {k: mp.mpf(v) for k,v in raw.items()}
t0=time.time()
li5v=[vals[f'li5:{i}'] for i in range(N)]
lg=[vals[f'log:{names[i]}'] for i in range(d)]
def combo_val(row, comp):
    v=sum(mp.mpf(c)*li5v[i] for i,c in enumerate(row) if c)
    for ms_s, cf_s in comp.items():
        ms=eval(ms_s); cf=Ff(cf_s)
        pv=mp.mpf(1)
        for ii in ms: pv*=lg[ii]
        v -= mp.mpf(cf.numerator)/cf.denominator*pv
    return v
combos=[combo_val(row,comp) for row,comp in zip(K,comps)]
keep=[j for j,v in enumerate(combos) if abs(v)>mp.mpf('1e-1300')]
print(f"combos: {len(combos)}, kept {len(keep)}, exact-zero {len(combos)-len(keep)} ({time.time()-t0:.0f}s)")
# zero-symbol product atoms
atoms=[]
for ms in itertools.combinations_with_replacement(range(d),3):
    v=vals['pi2']
    for i in ms: v*=lg[i]
    atoms.append((f"pi2*{ '*'.join(names[i] for i in ms) }", v))
for i in range(d):
    atoms.append((f"pi4*{names[i]}", vals['pi2']**2*lg[i]))
atoms.append(("z3*pi2", vals['z3']*vals['pi2']))
for a,b in itertools.combinations_with_replacement(range(d),2):
    atoms.append((f"z3*{names[a]}*{names[b]}", vals['z3']*lg[a]*lg[b]))
atoms.append(("T3*pi2", vals['T3']*vals['pi2']))
for a,b in itertools.combinations_with_replacement(range(d),2):
    atoms.append((f"T3*{names[a]}*{names[b]}", vals['T3']*lg[a]*lg[b]))
print(f"zero-symbol atoms: {len(atoms)}")
atom_names=['T','z5']+[f'K{j}' for j in keep]+[a[0] for a in atoms]
atom_vals=[vals['T'], vals['z5']]+[combos[j] for j in keep]+[a[1] for a in atoms]
print(f"LLL dimension: {len(atom_vals)}")
json.dump(dict(atom_names=atom_names, keep=keep), open('out/d13/atoms5.json','w'))
pickle.dump(dict(atom_names=atom_names, atom_vals=[mp.nstr(v,1440) for v in atom_vals]), open('out/d13/atoms5.pkl','wb'))
# timing calibration: LLL at 400 digits
import flint
def lll_at(P):
    mp.mp.dps=P+30
    n=len(atom_vals); Cc=mp.power(10,P)
    B=[]
    for i,v in enumerate(atom_vals):
        row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*Cc)); B.append(row)
    t1=time.time(); R=flint.fmpz_mat(B).lll(gram='exact'); dt=time.time()-t1
    # count short rows with a0!=0 and small residual
    mp.mp.dps=1450
    nrel=0; bestres=None
    for i in range(n):
        a=[int(R[i,j]) for j in range(n)]
        if not any(a) or a[0]==0: continue
        res=sum(mp.mpf(x)*v for x,v in zip(a,atom_vals))
        if abs(res)<mp.mpf(10)**(-(P-30)):
            nrel+=1
            if bestres is None or abs(res)<bestres: bestres=abs(res)
    print(f"LLL P={P}: dim {n}, {dt:.0f}s, relations(a0!=0)={nrel}, best residual {mp.nstr(bestres,3) if bestres else 'NA'}", flush=True)
    return nrel
lll_at(400)
