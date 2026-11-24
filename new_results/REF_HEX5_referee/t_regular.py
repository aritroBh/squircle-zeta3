import mpmath as mp, time, sys
from direct5 import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv)>1 else 30
r = 1/mp.sqrt(3)
V = [(r*mp.cos(mp.pi/2+i*mp.pi/3), r*mp.sin(mp.pi/2+i*mp.pi/3)) for i in range(3)]
H = setup(V, [1,1,1])
t0=time.time()
tot, per = all_derivs(H)
for k in range(6): print(k, mp.nstr(tot[k], mp.mp.dps))
print('levels', [p[2] for p in per], 'err', [mp.nstr(p[1],3) for p in per])
print('P(1)-2sqrt3', mp.nstr(tot[0]-2*mp.sqrt(3),3))
print('P5 - paper', mp.nstr(tot[5]-mp.mpf('656.4153619676247544766476090583112096482642646489854860418535509621431029'),3))
print('P3 - paper', mp.nstr(tot[3]-mp.mpf('-5.145455017135615837171681131123975071115553566419741721'),3))
print('time', time.time()-t0)
