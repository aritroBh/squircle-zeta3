import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import gengine, mpmath as mp
t0 = time.time()
G250, letters = gengine.gtab(250, 830)
print("dps-250 table built in %.1fs, %d words" % (time.time()-t0, len(G250)))
# quick validation at this precision
mp.mp.dps = 250
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
print("G(0^4,p;1)+Li_5(1/p):", mp.nstr(G250[(L0,)*4+(LP,)]+mp.polylog(5,1/LP), 4))
print("G(p,q;1) vs quad:", mp.nstr(G250[(LP,LQ)]-mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1]), 4))
