import sys, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc
import mpmath as mp
mp.mp.dps=680
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
d=pickle.load(open(WORK+'/cands29.pickle','rb'))
fw=d['freewords']; true29=d['cands']
t=time.time()
evsetup.setup(2500,4200,3800)   # prec 2500 bits (~752 d), N=3800 (~700 digits)
V1=[]
for w in fw:
    z=cc.h(tuple(w))
    V1.append(mp.mpc(mp.mpf(z.real.mid().str(700,radius=False,more=True)),
                     mp.mpf(z.imag.mid().str(700,radius=False,more=True))))
print(f'k=1 values at N=3800/prec2500: {time.time()-t:.0f}s')
def res(row):
    s=mp.mpc(0)
    for i,c in enumerate(row):
        if c: s+=c*V1[i]
    return abs(s)
print('--- 29 relations at ~700 digits ---')
worst=mp.mpf(0)
for i,row in enumerate(true29):
    r=res(row); worst=max(worst,r)
    print(f'{i:3d} log10res={float(mp.log10(r)) if r>0 else -999:8.1f}')
print('worst relation residual:', mp.nstr(worst,3))
# accident controls: reload the accident vectors from a fresh LLL? Instead recompute residuals
# for 5 accident combos: random integer combos of the 29 (should NOT vanish)
import random
random.seed(7)
print('--- random combos (controls) ---')
for trial in range(3):
    row=[sum(random.randint(-50,50)*true29[j][i] for j in range(29)) for i in range(110)]
    # this IS in the relation lattice, should vanish; use instead single words perturbed
    r=res(row)
    print(f'lattice-combo {trial}: log10res={float(mp.log10(r)) if r>0 else -999:8.1f} (should vanish)')
# non-relation control: single word
print('single-word control:', mp.nstr(res([1]+[0]*109),3),'(should be O(0.03-0.3))')
pickle.dump({'V1_700':[(mp.nstr(v.real,700),mp.nstr(v.imag,700)) for v in V1]},open(WORK+'/vals110_k1_700.pickle','wb'))
