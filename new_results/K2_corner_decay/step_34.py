
import mpmath as mp
from mpmath.calculus.quadrature import TanhSinh
mp.mp.dps = 100
alpha = mp.pi/3; c = mp.cos(alpha)
def B0(x):
    Z = mp.exp(-x)
    return (Z-2*c)/(1+mp.sqrt(1-2*c*Z+Z**2))
Tint0 = lambda x: B0(x)/(mp.exp(x)+1)
T0_quad = mp.quad(Tint0, [0, 1, 10, mp.inf])
print("T0 mp.quad =", mp.nstr(T0_quad, 30))
R0 = -mp.log(2)
M0 = mp.quad(lambda g: ((1+g)**(-1)*mp.sqrt(1+g*g-2*c*g)-1)/g, [0,1])
print("T0+R0-M0 (quad) =", mp.nstr(T0_quad+R0-M0, 6))
TS = TanhSinh(mp.mp)
for udeg in [5,6,7,8]:
    nodes=[]
    for d in range(1, udeg+1): nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
    uh = mp.mpf(2)**(-udeg)
    s = mp.fsum([w0*uh*Tint0(um) for um,w0 in nodes if um < 1200])
    print("udeg=%d nodes=%d  nodesum-T0quad = %.2e" % (udeg, len(nodes), s-T0_quad))
