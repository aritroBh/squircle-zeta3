import sys
sys.path.insert(0, '<project>/research/08_theory/SRF2_cube_residue')
sys.path.insert(0, '<project>/research/07_push/S3_superellipsoid/verify')
import mpmath as mp
from flint import arb, ctx
ctx.dps = 80; mp.mp.dps = 80
exec(open('<project>/research/08_theory/SRF2_cube_residue/step_37.py').read().split("print(\"J_reliable(-1)")[0])

print("J_reliable(-1) =", mp.nstr(J_reliable(arb(-1), mp.mpf(-1)), 25), "  target -0.84721308479397908661")
for eta in ['0.2','0.1','0.05','0.02','0.01','0.005','0.002']:
    e_mp = -1 + mp.mpf(eta)
    Jv = J_reliable(arb(mp.nstr(e_mp, 30)), e_mp)
    q = 1+2*e_mp
    print(f"eta={eta}: J={mp.nstr(Jv,22)}   -qJ = {mp.nstr(-q*Jv, 22)}")
