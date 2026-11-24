import time, sys, json
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=120
import mk_eval
for tag,a in [('2pi3',2*mp.pi/3),('pi3',mp.pi/3)]:
    t=time.time()
    M,nn=mk_eval.mk_all(a,80)
    json.dump([mp.nstr(x,50) for x in M], open(f'mk_{tag}_k80.json','w'))
    A=mp.sqrt(2*a*mp.sin(a))
    print(f"alpha={tag}: time={round(time.time()-t,1)}s nodes={nn} A={mp.nstr(A,10)}")
    for k in [40,60,80]:
        N=k+1; th=mp.sqrt(2*a*N); lead=-(A/N)*mp.e**(-th)*mp.cos(th)
        print(f"   M_{k}={mp.nstr(M[k],14)}  lead={mp.nstr(lead,14)}")
