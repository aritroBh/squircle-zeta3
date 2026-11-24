
import time
import gengine_dec
t0 = time.time()
G420, letters = gengine_dec.gtab(420, 1180, cachedir=".")
print("G-table dps420 N1180 built in %.1fs" % (time.time()-t0))
import mpmath as mp
mp.mp.dps = 420
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
phi = (1+mp.sqrt(5))/2
print("check G(0^4,-phi;1)+Li_5(-1/phi):", mp.nstr(abs(G420[(L0,)*4+(LP,)]+mp.polylog(5,1/LP)),3))
print("check G(0^4,phi^2;1)+Li_5(1/phi^2):", mp.nstr(abs(G420[(L0,)*4+(LQ,)]+mp.polylog(5,1/LQ**1)),3))
