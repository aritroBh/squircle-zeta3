import mpmath as mp, sys
from jet import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv)>1 else 50
U = [(mp.mpf(1),mp.mpf(0)),(mp.mpf(0),mp.mpf(1))]
for h in [mp.mpf(1)/32, mp.mpf(1)/64]:
    P = perimeter_jet(U, h=h)
    ex = [4*mp.sqrt(2), 0, mp.pi**2*mp.sqrt(2)/6, mp.sqrt(2)*(mp.zeta(3)-mp.pi**2/3)]
    print(h, [mp.nstr(P[k]-ex[k],5) for k in range(4)])
