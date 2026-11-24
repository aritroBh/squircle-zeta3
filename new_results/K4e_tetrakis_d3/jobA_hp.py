
import os, json, time
os.chdir("research/08_theory/K4e_tetrakis_d3")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 135
import k4c_driver
from k4c_driver import chamber_S
k4c_driver.WMIN = mpf('1e-95')
k4c_driver.LMIN = mpf('1e-85')
s2 = mp.sqrt(2)
NA = [[1/s2,1/s2,0],[1/s2,-1/s2,0],[1/s2,0,1/s2],[1/s2,0,-1/s2],[0,1/s2,1/s2],[0,1/s2,-1/s2]]
def log(*a):
    s=" ".join(str(x) for x in a)
    print(s, flush=True)
    open("bgAhp.log","a").write(s+"\n")
open("bgAhp.log","w").write("jobAhp start dps=%d\n"%mp.mp.dps)
S1 = 2*mp.sqrt(5)/3
for tag, HO, HI in [("h192hp", mpf(1)/192, mpf(1)/192), ("h256hp", mpf(1)/256, mpf(1)/256)]:
    t0=time.time()
    S, st = chamber_S(NA, 4, HO, HI, mpf('5.0'), mpf('5.0'), HI2=HI)
    el=time.time()-t0
    vals=[]
    for k in range(5):
        v=mp.re(S[k])
        vals.append(mp.nstr(v,140))
        log("[%s] [d^%d]S = %s"%(tag,k,mp.nstr(v,110)))
    log("[%s] el=%ds stats=%s"%(tag,el,st))
    log("[%s] d0 err vs 2sqrt5/3 = %s"%(tag,mp.nstr(abs(mp.re(S[0])-S1),3)))
    json.dump({"Sdk":vals,"stats":st,"HO":str(HO),"HI":str(HI)}, open("modelA_%s.json"%tag,"w"), indent=1)
log("jobAhp DONE")
