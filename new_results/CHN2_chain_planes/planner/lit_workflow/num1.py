import mpmath as mp, time
mp.mp.dps=30
u1,u2,u3,u4=[mp.mpf(x) for x in ('0.3','0.7','0.5','0.4')]
def G(u1,u2,u3): return mp.gamma(1+u2)*mp.gamma(1+u1+u3)/mp.gamma(1+u1+u2+u3)
t0=time.time()
V=G(u1,u2,u3)*mp.hyp3f2(-u4,u3,1+u1+u3,1+u3,1+u1+u2+u3,1)
print("V via 3F2",V,time.time()-t0)
mp.mp.dps=15
f=lambda s: s**(u2-1)*mp.quad(lambda t: t**(u3-1)*(1-s)**u1*(1-t)**u4,[0,1-s])
I=u2*u3*mp.quad(f,[0,1])
print("V via integral (15 dps)",I)
