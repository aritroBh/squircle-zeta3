import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
mp.mp.dps = 60
r2 = mp.sqrt(2); pval = 2+r2
t0=time.time()
# n=2: Li_2^real(p) = Re(polylog(2,p)) validated by single quad: -G(0,1/p;1) = int_0^1 log|1-tp|/t dt
v2 = mp.quad(lambda t: mp.log(abs(1-t*pval))/t, [0, 1/pval, 1], method='tanh-sinh')
print("step n=2: Re(polylog(2,p)) - int =", mp.nstr(mp.re(mp.polylog(2,pval)) - v2, 4))
# wait sign: Li_2^real(p) = -G(0,1/p;1) = -int_0^1 log|1-tp|/t dt? G(0,a;1) = int_0^1 log(1-t/a)/t dt with a=1/p: log(1-tp)
# = int_0^1 log|1-tp|/t dt (real). And Li_2^real(p) = -G(0,1/p;1) = -that.
print("   sign check: Re(polylog(2,p)) + int =", mp.nstr(mp.re(mp.polylog(2,pval)) + v2, 4))
# n=3: Li_3^real(p) = int_0^1 Li_2^real(tp)/t dt = int_0^1 Re(polylog(2,tp))/t dt
v3 = mp.quad(lambda t: mp.re(mp.polylog(2, t*pval))/t, [0, 1/pval, 1], method='tanh-sinh')
print("step n=3: Re(polylog(3,p)) - int =", mp.nstr(mp.re(mp.polylog(3,pval)) - v3, 4))
# n=4
v4 = mp.quad(lambda t: mp.re(mp.polylog(3, t*pval))/t, [0, 1/pval, 1], method='tanh-sinh')
print("step n=4: Re(polylog(4,p)) - int =", mp.nstr(mp.re(mp.polylog(4,pval)) - v4, 4))
# n=5
v5 = mp.quad(lambda t: mp.re(mp.polylog(4, t*pval))/t, [0, 1/pval, 1], method='tanh-sinh')
print("step n=5: Re(polylog(5,p)) - int =", mp.nstr(mp.re(mp.polylog(5,pval)) - v5, 4))
print("total %.1fs" % (time.time()-t0))
