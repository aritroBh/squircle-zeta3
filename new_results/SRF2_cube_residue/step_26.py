import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
from vlib import ts_nodes
from srf2_lib import W_int, J_int, tomp
ctx.dps = 60; mp.mp.dps = 60

eps = arb('0.1')
q = 1+2*eps
nodes = ts_nodes(arb(1)/32, 4.6)

# direct raw V(0.1) = int int a^{2e-1} c^{e-1} (g - g0)
tot = arb(0)
for (a, am, wa, pa) in nodes:
    la = a.log(); l1a = (1+a).log()
    A2 = a*a*(arb(-2*eps)*la).exp()
    S0 = (1+A2).sqrt()
    g0 = (arb(-2*eps)*l1a).exp()/(1+a)*S0
    prefa = (arb(2*eps)*la).exp()/a
    row = arb(0)
    for (c, cm, wc, pc) in nodes:
        lc = c.log()
        C2 = c*c*(arb(-2*eps)*lc).exp()
        N = 1+a+a*c
        g = (arb(-2*eps)*N.log()).exp()/N*(1+A2+A2*C2).sqrt()
        prefc = (arb(eps)*lc).exp()/c
        row += prefc*(g-g0)*wc
    tot += prefa*row*wa
print("direct V(0.1) =", mp.nstr(tomp(tot), 20))

# direct int int w * c B1
tot2 = arb(0)
for (a, am, wa, pa) in nodes:
    la = a.log(); l1a = (1+a).log()
    A2 = a*a*(arb(-2*eps)*la).exp()
    S0 = (1+A2).sqrt()
    B1 = -q*a*(arb(-2*eps)*l1a).exp()/(1+a)**2*S0
    prefa = (arb(2*eps)*la).exp()/a
    row = arb(0)
    for (c, cm, wc, pc) in nodes:
        lc = c.log()
        prefc = (arb(eps)*lc).exp()/c
        row += prefc*c*wc
    tot2 += prefa*B1*row*wa
print("direct int w cB1 =", mp.nstr(tomp(tot2), 20))

W = W_int(eps); J = J_int(eps)
print("W(0.1) =", mp.nstr(tomp(W), 20))
print("J(0.1) =", mp.nstr(tomp(J), 20))
print("-qJ/(e+1) =", mp.nstr(tomp(-q*J/(eps+1)), 20))
print("W + direct cB1 =", mp.nstr(tomp(W + tot2), 20))
