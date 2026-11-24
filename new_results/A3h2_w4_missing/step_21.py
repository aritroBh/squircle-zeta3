import mpmath as mp
mp.mp.dps=120
def a2mp(z, digits=100):
    return mp.mpc(mp.mpf(z.real.mid().str(digits,radius=False,more=True)),
                  mp.mpf(z.imag.mid().str(digits,radius=False,more=True)))
z3=cc.h((-1,-1,0))
print('zeta3 err:', mp.nstr(abs(a2mp(z3)-mp.zeta(3)),5))
# check convergence: N=800 vs N=400
cc.N=400
for f in (cc.series,cc.half,cc.g_at_one,cc.h): f.cache_clear()
z3b=cc.h((-1,-1,0))
z2b=cc.h((7,4,7,7))
print('zeta3 N400 err:', mp.nstr(abs(a2mp(z3b)-mp.zeta(3)),5))
cc.N=800
for f in (cc.series,cc.half,cc.g_at_one,cc.h): f.cache_clear()
z3a=cc.h((-1,-1,0)); z2a=cc.h((7,4,7,7))
print('N800 vs N400 zeta3:', mp.nstr(abs(a2mp(z3a)-a2mp(z3b)),5))
print('N800 vs N400 (7,4,7,7):', mp.nstr(abs(a2mp(z2a)-a2mp(z2b)),5))
