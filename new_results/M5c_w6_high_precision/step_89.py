import sys, pickle, json, math
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3620
from fractions import Fraction as F
from w4_symbol import Q, R, A, B
from m5b_common import pool_zs
from w6_solve import atoms_w6, BASE_LI6

pool=pool_zs()
li6=list(BASE_LI6)+[pool[i] for i in [14,15,13,9,7,6]]
basis=atoms_w6(li6)
sol=json.load(open(work+'/w6x14_15_13_9_7_6_symbol_solve.json'))
K=sol['Ksat']
vals=[mp.mpf(v) for v in pickle.loads(open(work+'/w6_atomvals_3500_x14_15_13_9_7_6.pkl','rb').read())]
assert len(vals)==len(basis)==155
kv=[mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
# which kernel vectors involve the extra Li6 atoms (indices 149..154)?
extra_names=[b[0] for b in basis[149:]]
print('extra atoms:', extra_names)
small=[]
for j,krow in enumerate(K):
    nz=[(basis[i][0],krow[i]) for i in range(155) if krow[i] and i>=149]
    if nz:
        print(f'kv[{j}] involves {nz}, |value|={mp.nstr(abs(kv[j]),4)}')
        if abs(kv[j])<mp.mpf(10)**(-3000): small.append(j)
print('exactly-zero kv columns (|v|<1e-3000):', small)
