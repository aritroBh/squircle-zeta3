
import os, json, time
os.chdir("research/08_theory/K4e_tetrakis_d3")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 160
import k4c_driver
from k4c_driver import chamber_S
k4c_driver.WMIN = mpf('1e-125')
k4c_driver.LMIN = mpf('1e-115')
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
def log(*a):
    s=" ".join(str(x) for x in a)
    print(s, flush=True)
    open("bgA384.log","a").write(s+"\n")
open("bgA384.log","w").write("jobA384 start dps=%d\n"%mp.mp.dps)
S1 = 2*mp.sqrt(5)/3
for tag, HO, HI in [("h384", mpf(1)/384, mpf(1)/384)]:
    t0=time.time()
    S, st = chamber_S(NA, 4, HO, HI, mpf('5.4'), mpf('5.4'), HI2=HI)
    el=time.time()-t0
    vals=[]
    for k in range(5):
        v=mp.re(S[k])
        vals.append(mp.nstr(v,160))
        log("[%s] [d^%d]S = %s"%(tag,k,mp.nstr(v,150)))
    log("[%s] el=%ds stats=%s"%(tag,el,st))
    log("[%s] d0 err vs 2sqrt5/3 = %s"%(tag,mp.nstr(abs(mp.re(S[0])-S1),3)))
    json.dump({"Sdk":vals,"stats":st,"HO":str(HO),"HI":str(HI)}, open("modelA_%s.json"%tag,"w"), indent=1)
log("jobA384 DONE")
