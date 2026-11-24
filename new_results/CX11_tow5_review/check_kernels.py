"""Independent exact full-kernel Newton checks; no submitted verifier imported."""
import pickle
from fractions import Fraction as F
from math import comb
from pathlib import Path
root=Path('../WTS_weight_k_structure/theory')
K={k:pickle.load(open(root/f'kern_{k}.pkl','rb')) for k in range(6,14)}
keys=[(6,0),(5,1),(4,2),(3,3)]; z=(F(0),)*4
sparse=bad6=0; full=0; coefficient_kernels={}
for key in keys:
 poles=set().union(*(set(K[k][key]['real']) for k in K))
 assert poles <= {(F(0),F(0)),(F(1),F(1)),(F(-1),F(-1)),(F(1),F(-1))}
 # Use k=8..13 finite differences and evaluate the Newton polynomial at every k.
 coef={}
 for pole in sorted(poles):
  coef[pole]=[]
  for c in range(4):
   col=[K[k][key]['real'].get(pole,z)[c] for k in range(8,14)]
   ds=[]
   while col:
    ds.append(col[0]); col=[b-a for a,b in zip(col,col[1:])]
   def choose(n,j):
    a=F(1)
    for h in range(j): a*=F(n-h,h+1)
    return a
   def pred(k): return sum(ds[j]*choose(k-8,j) for j in range(6))
   for k in range(7,14):
    assert pred(k)==K[k][key]['real'].get(pole,z)[c],(key,pole,c,k)
    full+=1
   bad6+=pred(6)!=K[6][key]['real'].get(pole,z)[c]
   sparse+=1
   coef[pole].append(ds)
 for k in range(6,14):
  item=K[k][key]
  assert set(item)=={'real','ipos','ineg','const'}
  assert not item['ipos'] and not item['ineg'] and item['const']==z,(key,k)
 # Produce the six coefficient kernels based at k=7; integration defines E_{6,m}.
 coefficient_kernels[key]={}
 for pole in sorted(poles):
  column=[K[k][key]['real'].get(pole,z) for k in range(7,13)]
  cs=[]
  while column:
   cs.append(column[0]);column=[tuple(b[c]-a[c] for c in range(4)) for a,b in zip(column,column[1:])]
  coefficient_kernels[key][pole]=cs
print(f'All four keys: {sparse} sparse real QSI coordinates, {full} exact evaluations k=7..13 PASS.')
print(f'k=6: {bad6}/{sparse} coordinates off the same quintic; all +/-i and constant parts zero at k=6..13.')
assert sparse==56 and bad6==18
pickle.dump(coefficient_kernels,open('W7_Newton_kernels.pkl','wb'))
print('Six forward-difference coefficient kernels based at k=7 saved. Fixed integration gives E_{6,1..6}.')
