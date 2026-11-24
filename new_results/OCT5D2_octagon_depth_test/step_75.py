import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import gengine, mpmath as mp
t0 = time.time()
G430, letters = gengine.gtab(430, 1440)
print("dps-430 table built in %.1fs" % (time.time()-t0))
mp.mp.dps = 430
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
print("G(0^4,p;1)+Li_5(1/p):", mp.nstr(G430[(L0,)*4+(LP,)]+mp.polylog(5,1/LP), 4))
print("G(p,q;1) vs quad:", mp.nstr(G430[(LP,LQ)]-mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1]), 4))
