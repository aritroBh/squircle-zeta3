import os, json, time
os.chdir("research/08_theory/K4d_identify")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 90
from k4c_driver import chamber_S
s3 = 1/mp.sqrt(3)
NB = [[1,0,0],[0,1,0],[0,0,1],[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
def log(*a):
    s=" ".join(str(x) for x in a)
    print(s, flush=True)
    open("bgB.log","a").write(s+"\n")
open("bgB.log","w").write("jobB start dps=%d\n"%mp.mp.dps)
import sympy as sp
S1 = mpf(sp.N(-36 + 21*sp.sqrt(3) + sp.Rational(9,13)*sp.sqrt(35515-20504*sp.sqrt(3)), 95))
log("exact S1 = "+mp.nstr(S1,70))
for tag, HO, HI, HI2 in [("h48", mpf(1)/48, mpf(1)/64, mpf(1)/64), ("h96o", mpf(1)/96, mpf(1)/64, mpf(1)/64)]:
    t0=time.time()
    S, st = chamber_S(NB, 4, HO, HI, mpf('4.5'), mpf('4.5'), HI2=HI2)
    el=time.time()-t0
    vals=[]
    for k in range(5):
        v=S[k]
        im = abs(mp.im(v)) if mp.im(v)!=0 else mpf(0)
        v = mp.re(v)
        vals.append(mp.nstr(v,100))
        log("[%s] [d^%d]S = %s  |imag|=%s"%(tag,k,mp.nstr(v,70),mp.nstr(im,2)))
    log("[%s] el=%ds stats=%s"%(tag,el,st))
    log("[%s] d0 err vs exact = %s"%(tag,mp.nstr(abs(mp.re(S[0])-S1),3)))
    json.dump({"Sdk":vals,"stats":st,"HO":str(HO),"HI":str(HI),"HI2":str(HI2)}, open("modelB_%s.json"%tag,"w"), indent=1)
log("jobB DONE")
