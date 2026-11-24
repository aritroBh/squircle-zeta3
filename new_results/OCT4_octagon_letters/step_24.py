import os, time, sys; os.chdir('research/08_theory/OCT4_octagon_letters')
sys.path.insert(0,'research/08_theory/D15b_polygon_law')
import mpmath as mp
import poly_coeffs as pc
from oct4_common import log
# gauge timing at moderate dps
for dps in [80, 120]:
    mp.mp.dps = dps
    t0=time.time()
    r = pc.edge_arc_coeffs(4, kmax=3)
    dt=time.time()-t0
    print(f"dps={dps}: P'''(1) = {mp.nstr(r[3],30)}  ({dt:.1f}s)")
mp.mp.dps=120
log("step: quadrature timing test done; P'''(1) framework OK")
