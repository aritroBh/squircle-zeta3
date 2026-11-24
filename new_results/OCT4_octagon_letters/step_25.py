import os, time, sys; os.chdir('research/08_theory/OCT4_octagon_letters')
sys.path.insert(0,'research/08_theory/D15b_polygon_law')
import mpmath as mp
import poly_coeffs as pc
from oct4_common import log
mp.mp.dps = 300
t0=time.time()
D = pc.edge_arc_coeffs(4, kmax=3)
P3 = 6*D[3]
u = mp.sqrt(2+mp.sqrt(2))
print("P'''(1) =", mp.nstr(P3, 60))
print("time:", round(time.time()-t0,1),"s")
print("P'''/u =", mp.nstr(P3/u, 60))
# save high precision targets
import json
targets = {'dps': mp.mp.dps, 'P3': mp.nstr(P3, mp.mp.dps-5), 'P3_over_u': mp.nstr(P3/u, mp.mp.dps-5),
           'P0': mp.nstr(D[0],60), 'P1': mp.nstr(D[1],60), 'P2': mp.nstr(2*D[2],60)}
json.dump(targets, open('targets300.json','w'))
log("step: P'''(1) and P'''/u at %d dps; P'''=%s ; saved targets300.json"%(mp.mp.dps, mp.nstr(P3,25)))
