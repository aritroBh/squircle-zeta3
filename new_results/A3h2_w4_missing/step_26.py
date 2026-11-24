import sys, json, pickle
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, h_k, sig, a2mp
from flint import acb
import mpmath as mp
mp.mp.dps=150
evsetup.setup(900,1200,1000)
BASE='<project>'
rogers=json.load(open(BASE+'/research/08_theory/A3f_level8_relations/weight2_new_relation.json'))['row']
def evalrel_acb(terms,k):
    tot=acb(0)
    for w,c in terms: tot += int(c)*h_k(w,k)
    return tot
v3=evalrel_acb([(tuple(w),c) for w,c in rogers],3)
z=a2mp(v3,140)
print('five-term k=3 residual:', mp.nstr(z,40))
print('imag/pi:', mp.nstr(z.imag/mp.pi,20))
print('imag/(pi*sqrt2):', mp.nstr(z.imag/(mp.pi*mp.sqrt(2)),20))
print('imag/(2*pi*log2):', mp.nstr(z.imag/(2*mp.pi*mp.log(2)),20))
print('imag/(pi*log2):', mp.nstr(z.imag/(mp.pi*mp.log(2)),20))
print('imag/(pi**2):', mp.nstr(z.imag/mp.pi**2,20))
# PSLQ to identify: candidates pi*ln2, pi^2, pi*sqrt2...
consts=[mp.pi*mp.log(2), mp.pi**2, mp.pi*mp.sqrt(2), mp.pi*mp.log(1+mp.sqrt(2)), mp.pi]
import itertools
best=None
for c in consts:
    q=z.imag/c
    print(mp.nstr(c,12),'->',mp.nstr(q,18))
