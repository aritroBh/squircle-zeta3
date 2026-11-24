import sys, pickle
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, a2mp
import mpmath as mp
mp.mp.dps=160
BASE='<project>'
vh=pickle.load(open(BASE+'/research/08_theory/A3f4_glanois_basis/vals_hi.pickle','rb'))
print(vh.keys(), vh['settings'])
print('n vals:', len(vh['vals']), 'freewords:', len(vh['freewords']))
# my evaluation of the same 32 free words
evsetup.setup(1200,1500,1000)
worst=mp.mpf(0); worstw=None
for i,w in enumerate(vh['freewords']):
    z=a2mp(cc.h(tuple(w)),150)
    re_s,im_s=vh['vals'][27+i]
    mp.mp.dps=160
    zref=mp.mpc(mp.mpf(re_s),mp.mpf(im_s))
    err=abs(z-zref)
    if err>worst: worst=err; worstw=w
print('max |mine - A3f4 vals_hi| over 32 words:', mp.nstr(worst,4), 'at', worstw)
