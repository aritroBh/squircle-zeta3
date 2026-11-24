import sys; sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup, mpmath as mp
mp.mp.dps=180
# Convergence certificate at prec 700 (210 digits): N=800 vs N=400
evsetup.setup(700,1000,800); za=evsetup.cc.h((7,4,7,7)); z3a=evsetup.cc.h((-1,-1,0))
evsetup.setup(700,1000,400); zb=evsetup.cc.h((7,4,7,7)); z3b=evsetup.cc.h((-1,-1,0))
print('|N800-N400| (7,4,7,7):', mp.nstr(abs(evsetup.a2mp(za,170)-evsetup.a2mp(zb,170)),5))
print('|N800-N400| zeta3   :', mp.nstr(abs(evsetup.a2mp(z3a,170)-evsetup.a2mp(z3b,170)),5))
# internal accuracy of values: radius
print('radius (7,4,7,7):', za.real.rad().str(3,radius=False), za.imag.rad().str(3,radius=False))
# worst-letter probe: (7,7,7,7) is depth4 all 7; (5,) ends...
for w in [(7,7,7,7),(7,6,3,3)]:
    evsetup.setup(700,1000,800); xa=evsetup.cc.h(w)
    evsetup.setup(700,1000,400); xb=evsetup.cc.h(w)
    print(w,'|N800-N400|:', mp.nstr(abs(evsetup.a2mp(xa,170)-evsetup.a2mp(xb,170)),5))
