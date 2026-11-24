
import time
import gengine_dec
t0 = time.time()
G250, letters = gengine_dec.gtab(250, 700, cachedir=".")
print("G-table dps250 N700 built in %.1fs, %d words" % (time.time()-t0, len(G250)))
import mpmath as mp
mp.mp.dps = 250
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
phi = (1+mp.sqrt(5))/2
print("check G(0^4,-phi;1)+Li_5(-1/phi):", mp.nstr(abs(G250[(L0,)*4+(LP,)]+mp.polylog(5,1/LP)),3))
