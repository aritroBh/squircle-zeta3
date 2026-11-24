import mpmath as mp, time
from direct import *
mp.mp.dps = 30
r3 = 1/mp.sqrt(3)
V = [(r3*mp.cos(mp.pi/6), r3*mp.sin(mp.pi/6)), (0, r3), (r3*mp.cos(5*mp.pi/6), r3*mp.sin(5*mp.pi/6))]
H = setup(V, [1,1,1])
print('c', H['c'])
t=time.time()
tot=[0]*4
for q in range(3):
    e = edge_jet_fast(H, q)
    tot=[tot[i]+2*e[i] for i in range(4)]
print(time.time()-t)
print('P(1)', tot[0], 2*mp.sqrt(3))
print("P'(1)", tot[1], 2*mp.sqrt(3)*mp.log(2))
l=mp.log(2)
print("P'''", tot[3], mp.sqrt(3)/8*(16*l**3+3*mp.pi**2*l-96*l**2+96*l-9*mp.zeta(3)-6*mp.pi**2))
