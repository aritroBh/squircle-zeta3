import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp
mp.mp.dps = 60
# (a) J = int_0^inf (s g1^2/2 + g1^3/6) ds =? zeta(4)/32
g1 = lambda s: mp.log(1+mp.exp(-2*s))
J1 = mp.quad(lambda s: s*g1(s)**2/2, [0, mp.inf])
J2 = mp.quad(lambda s: g1(s)**3/6, [0, mp.inf])
J = J1 + J2
print("J =", mp.nstr(J, 50))
print("zeta(4)/32 =", mp.nstr(mp.zeta(4)/32, 50), " resid:", mp.nstr(J-mp.zeta(4)/32, 6))
# (b) F1 A4 PSLQ on weight-4 constants
A4 = mp.mpf("0.27059026108183716150330571099319")
li4h = mp.polylog(4, mp.mpf('0.5'))
basis = [mp.zeta(4), li4h, mp.log(2)**4, mp.pi**2*mp.log(2)**2, mp.zeta(4)*mp.log(2), mp.pi**2, mp.log(2)**2, mp.log(2), mp.mpf(1)]
M = [A4] + basis
rel = mp.pslq(M, tol=mp.mpf('1e-30'), maxcoeff=10**10, maxsteps=10000)
print("F1 A4 PSLQ:", rel)
if rel: print("  resid:", mp.nstr(mp.fsum(mp.mpf(r)*x for r,x in zip(rel,M)),6))
