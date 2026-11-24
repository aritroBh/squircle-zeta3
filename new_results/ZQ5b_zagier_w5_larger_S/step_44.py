import os, json
from fractions import Fraction
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
step = r'''
import sys, json, time
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
ZQ='../ZQ5_zagier_weight5'
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),17), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
S = json.load(open('d17_min_setup.json'))
names=S['names']; N=S['N']; M=S['M']
prod_list=[(p[0],tuple(p[1]),p[2]) for p in S['prod_list']]
P=1200
mp.mp.dps=P+30
t0=time.time()
s=mp.sqrt(17)
lg={n: mp.log(abs(g.emb(P+30))) for n,g in zip(names, W.make_field(17,split=[2],inert=[])['letters'])}
vals={}
vals['z5']=mp.zeta(5)
vals['T']=W.Lchi_hurwitz(17,5,P)/s
vals['pi2']=mp.pi**2
for i,c in enumerate(cands):
    vals[f'li5:{i}']=W.li5(c['x'].emb(P+30),P+10)
for kind,ms,nm in prod_list:
    if kind=="log5":
        v=mp.mpf(1)
        for i in ms: v*=lg[names[i]]
    elif kind=="pi2log3":
        v=vals['pi2']
        for i in ms: v*=lg[names[i]]
    else:
        v=vals['pi2']**2*lg[names[ms[0]]]
    vals[f'prod:{nm}']=v
# store
store = {k: mp.nstr(v, P+25) for k,v in vals.items()}
import pickle
pickle.dump(store, open('d17_vals_1200.pkl','wb'))
print(f"values computed at {P} digits ({time.time()-t0:.0f}s)")
# verify relation residual
mp.mp.dps=P+30
c_rel=[Ff(x) for x in S['c_rel']]
tot = Ff(S['c_rel'][0])*0  # placeholder
res = c_rel[0]*vals['z5']
for i in range(N):
    if c_rel[1+i]: res += c_rel[1+i]*vals[f'li5:{i}']
for pi in range(M):
    if c_rel[1+N+pi]: res += c_rel[1+N+pi]*vals[f'prod:{prod_list[pi][2]}']
# c_rel is coefficient of RHS terms; relation is T_coef*T = sum -> residual = T_coef*T - RHS
Tres = Ff(-1419857)*vals['T'] - res
print("relation residual (T_coef*T - RHS):", mp.nstr(Tres,4))
# verify zero rows
worst=mp.mpf(0)
for zi,row in enumerate(S['Z']):
    z=[Ff(x) for x in row]
    r = z[0]*vals['z5']
    for i in range(N):
        if z[1+i]: r += z[1+i]*vals[f'li5:{i}']
    for pi in range(M):
        if z[1+N+pi]: r += z[1+N+pi]*vals[f'prod:{prod_list[pi][2]}']
    if abs(r)>worst: worst=abs(r); wzi=zi
print(f"worst zero-row residual over {len(S['Z'])} rows: {mp.nstr(worst,4)} (row {wzi})")
'''
open('step_04b_d17_vals.py','w').write(step)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_04b_d17_vals.py'], capture_output=True, text=True, cwd='.', timeout=2400)
print(r.stdout); print("ERR:", r.stderr[-1500:])
