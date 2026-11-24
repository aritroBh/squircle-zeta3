
from flint import ctx, arb
ctx.dps = 40
from dcy2_cell import *
lo, hi = 83966.56115454846, 85645.89237763942
al = alpha_of(lo); c_lo, c_hi = c2min_of(lo), c2min_of(hi)
D_lo, D_hi = Delta_of(lo), Delta_of(hi)
gam = al
gamD_lo = gam*D_lo
print("gamD_lo =", gamD_lo.str(10), " c_lo/2 =", (c_lo/2).str(10))
if gamD_lo <= c_lo/2:
    print("branch 1")
else:
    print("branch 2")
    T = min(D_lo, c_lo/(2*gam))
    gamT3 = gam*T**3
    print("T =", T.str(8), " gamT3 =", gamT3.str(8))
    base = gamT3.exp()*Mtilde_of(hi, 12)
    print("base =", base.str(8))
    print("alpha^4/24*base =", (al**4/24*base).str(8))
    excess_hi = gamD_lo - c_hi
    print("excess_hi =", excess_hi.str(8))
    e2 = (ReL1max_of(lo)*D_lo).exp()
    if excess_hi < 0:
        rate_lo = c_hi - gamD_lo
        edge = D_hi**(12+1)*(-rate_lo*T**2).exp()*e2/(12+1)
    else:
        edge = e2*(excess_hi*D_hi**2).exp()*D_hi**(12+1)/(12+1)
    print("edge =", edge.str(8))
    print("GsplitB total =", (base+edge).str(8))
    print("alpha^4/24*(base+edge) =", (al**4/24*(base+edge)).str(8))
