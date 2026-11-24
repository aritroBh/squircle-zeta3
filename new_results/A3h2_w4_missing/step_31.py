import sys, json, time, pickle
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, sig, a2mp
import mpmath as mp
mp.mp.dps=380
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
freewords=out['freewords']
# per-word certificate N=1400 vs stored N=2000
worst={}
for k in (1,3,5,7):
    evsetup.setup(1200,2400,1400)
    wk=mp.mpf(0); ww=None
    for w in freewords:
        lo=cc.h(sig(w,k))
        re_s,im_s=out['vals'][k][w]
        hi=mp.mpc(mp.mpf(re_s),mp.mpf(im_s))
        d=abs(a2mp(lo,340)-hi)
        if d>wk: wk=d; ww=w
    worst[k]=(wk,ww)
    print(f'k={k}: max|N2000-N1400|={mp.nstr(wk,4)} at {ww}')
out['cert1400']={k:mp.nstr(v[0],4) for k,v in worst.items()}
pickle.dump(out,open(WORK+'/vals110.pickle','wb'))
