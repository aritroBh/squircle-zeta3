
import sys, json, pickle
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
S = json.load(open('d17_min_setup.json'))
N=S['N']; M=S['M']
prod_list=[(p[0],tuple(p[1]),p[2]) for p in S['prod_list']]
raw = pickle.load(open('d17_vals_1200.pkl','rb'))
mp.mp.dps = 1250
vals = {k: mp.mpf(v) for k,v in raw.items()}
c_rel=[Ff(x) for x in S['c_rel']]
res = c_rel[0]*vals['z5']
for i in range(N):
    if c_rel[1+i]: res += c_rel[1+i]*vals[f'li5:{i}']
for pi in range(M):
    if c_rel[1+N+pi]: res += c_rel[1+N+pi]*vals[f'prod:{prod_list[pi][2]}']
Tres = Ff(-1419857)*vals['T'] + res
print("relation residual a1[0]*T + RHS:", mp.nstr(Tres,4))
