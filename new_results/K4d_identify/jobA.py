import os, json, time
os.chdir("research/08_theory/K4d_identify")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
from k4c_driver import chamber_S
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
def log(*a):
    s=" ".join(str(x) for x in a)
    print(s, flush=True)
    open("bgA.log","a").write(s+"\n")
open("bgA.log","w").write("jobA start dps=%d\n"%mp.mp.dps)
S1 = 2*mp.sqrt(5)/3
for tag, HO, HI in [("h128", mpf(1)/128, mpf(1)/128), ("h192", mpf(1)/192, mpf(1)/192)]:
    t0=time.time()
    S, st = chamber_S(NA, 4, HO, HI, mpf('4.5'), mpf('4.5'), HI2=HI)
    el=time.time()-t0
    vals=[]
    for k in range(5):
        v=S[k]
        im = abs(mp.im(v)) if mp.im(v)!=0 else mpf(0)
        v = mp.re(v)
        vals.append(mp.nstr(v,100))
        log("[%s] [d^%d]S = %s  |imag|=%s"%(tag,k,mp.nstr(v,70),mp.nstr(im,2)))
    log("[%s] el=%ds stats=%s"%(tag,el,st))
    log("[%s] d0 err vs 2sqrt5/3 = %s"%(tag,mp.nstr(abs(mp.re(S[0])-S1),3)))
    json.dump({"Sdk":vals,"stats":st,"HO":str(HO),"HI":str(HI)}, open("modelA_%s.json"%tag,"w"), indent=1)
log("jobA DONE")
