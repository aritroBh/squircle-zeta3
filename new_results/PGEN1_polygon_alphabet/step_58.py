import time, pickle
open("progress.log","a").write(f"""{time.strftime('%Y-%m-%d %H:%M:%S')} steps 10-32:
- machinery (letters, wedge, delta) + ALL CONTROLS PASS (delta(1/2)=0, Landen triple delta=0 on 7 samples, hexagon args in CR for 6 r's, delta=0 at r=2,1/2)
- generic octagon: p_j per edge (e0:-2/3,12/7; e1:-3,-9; e2:-3/5,9/5; e3:-4,-24); CR_e (30,30,30,18 pts)
- explicit CR-only Li3 presentations of all four TW3_e (PROVED exact symbol identities, diff in D)
- delta(N3) per edge from OCT9d basis == from CR basis (PROVED identical); NONZERO for every edge
- regular octagon delta != 0: (r2^u)(x)r2(-5/2+15sqrt2/8) + (r2^u)(x)u(-11/4+2sqrt2)
- decagon delta = 0 (letters collapse to span(phi)); canonical pair pending
- 12-gon delta != 0 (9+9 terms over letters g2,g3,u); invariant under duplication presentation change (z4=z1^2)
- wedge3(delta)=0 (Lie^3) for octagon and 12-gon
- alphabet checks: regular octagon 6 pts in CR (2 orbits + orbit(1/2)); decagon CR = single orbit containing +-1/phi; 12-gon |CR|=63: 7 of 8 DOD5 pts in CR, z4=z1^2 duplication; span s(CR12) mod D = full 8-dim; all 8 symbols in span
""")
import flint
print("flint version:", flint.__version__)
from flint import fmpq_mat, fmpz_mat
print("fmpq_mat OK")
