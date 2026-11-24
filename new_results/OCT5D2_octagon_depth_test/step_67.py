import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 60
r2 = mp.sqrt(2); pval = 2+r2
# Li_4^real(p) = -G(0,0,0,1/p;1) = -int_0^1 dt/t int_0^t ds/s int_0^s dr/r log|1-r p|
nested = mp.quad(lambda t: mp.quad(lambda s: mp.quad(lambda r: mp.log(abs(1-r*pval))/r, [0,s])/(s), [0,t])/(t), [0,1], method='tanh-sinh')
print("Re(polylog(4,p)) + nested =", mp.nstr(mp.re(mp.polylog(4,pval))+nested, 4))
# same for Li_5: -G(0,0,0,0,1/p;1)
nested5 = mp.quad(lambda t: mp.quad(lambda s: mp.quad(lambda r: mp.quad(lambda u: mp.log(abs(1-u*pval))/u, [0,r])/r, [0,s])/(s), [0,t])/(t), [0,1], method='tanh-sinh')
print("Re(polylog(5,p)) + nested5 =", mp.nstr(mp.re(mp.polylog(5,pval))+nested5, 4))
