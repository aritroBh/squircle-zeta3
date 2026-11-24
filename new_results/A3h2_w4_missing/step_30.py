import sys, json, time, pickle
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, sig, a2mp
import mpmath as mp
mp.mp.dps=400
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
freewords=[tuple(w) for w in json.load(open(BASE+'/research/08_theory/A3h_weight4_target/T3_quotient_nf.json'))['freewords']]
print('n free words:', len(freewords))
def s350(z):
    return (z.real.mid().str(350,radius=False,more=True), z.imag.mid().str(350,radius=False,more=True))
# hi run: prec=1200 bits (~360 d), N=2000 ; lo: N=1000
out={'freewords':freewords,'settings':{'prec':1200,'cap':2400,'N_hi':2000,'N_lo':1000},'vals':{},'cert':{}}
t0=time.time()
for k in (1,3,5,7):
    evsetup.setup(1200,2400,2000)
    hi={}
    for i,w in enumerate(freewords):
        hi[w]=cc.h(sig(w,k))
    evsetup.setup(1200,2400,1000)
    worst=mp.mpf(0); worstw=None
    for w in freewords:
        lo=cc.h(sig(w,k))
        d=abs(a2mp(hi[w],340)-a2mp(lo,340))
        if d>worst: worst=d; worstw=w
    out['vals'][k]={w:s350(hi[w]) for w in freewords}
    out['cert'][k]=mp.nstr(worst,4)
    print(f'k={k}: max|N2000-N1000|={mp.nstr(worst,4)} at {worstw}; elapsed={time.time()-t0:.0f}s',flush=True)
pickle.dump(out,open(WORK+'/vals110.pickle','wb'))
print('saved vals110.pickle; total',time.time()-t0,'s')
