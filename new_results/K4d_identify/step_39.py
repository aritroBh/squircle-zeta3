import mpmath as mp, json, os, time
from mpmath import mpf
mp.mp.dps=70
WORK='research/08_theory/K4d_identify'
print("files:", [f for f in os.listdir(WORK) if f.endswith('.json')])
def digs(a,b):
    d=abs(a-b)
    return float(-mp.log10(d)) if d>0 else 999.0
SA96=[mpf(s) for s in json.load(open(WORK+'/modelA_h96.json'))['Sdk']]
SA128=[mpf(s) for s in json.load(open(WORK+'/modelA_h128.json'))['Sdk']]
print("\nModel A h96 vs h128 digit agreement:")
for k in range(5):
    print("  d^%d: %.1f digits"%(k, digs(SA96[k],SA128[k])))

# quick PSLQ speed test: recover q_A = 13/15 + (1/5)ln(3^13/2^6) i.e. 15 q_A = 13 -18 ln2 +39 ln3
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
qA = SA128[1]/SA128[0]
t0=time.time()
rel,resid = pslq([qA, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-50'), maxsteps=300)
print("\nPSLQ speed test (4-dim): time=%.2fs"%(time.time()-t0))
print("  qA relation:", rel, " resid:", mp.nstr(resid,3))
print("  expect ~ [15,-13,18,-39] up to scale")
