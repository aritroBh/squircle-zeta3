
import mpmath as mp
from mpmath.calculus.quadrature import TanhSinh
mp.mp.dps = 60
alpha = mp.pi/2; c = mp.cos(alpha)

def G0(x):  # T-node + R-node at q=0
    ex = mp.exp(-x)
    EZ = ex
    SQ = mp.sqrt(1-2*c*EZ+EZ**2)
    Bt0 = (1/ex)*(SQ-1)
    return Bt0/(mp.exp(x)+1) + (1/(1+ex)-1)

v_quad = mp.quad(G0, [0,1,10,mp.inf])
print("mp.quad int G0 =", mp.nstr(v_quad, 25), "  (should be M0=-0.553303...)")

TS = TanhSinh(mp.mp)
for udeg in [5,6]:
    nodes=[]
    for d in range(1, udeg+1): nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
    uh = mp.mpf(2)**(-udeg)
    tot = mp.mpf(0)
    for um, w0 in nodes:
        w = w0/(4*mp.cosh(um/2)**2)
        tot += w*uh*G0(um)
    print("udeg=%d  node-sum = %s   err=%.2e" % (udeg, mp.nstr(tot,25), abs(tot-v_quad)))
