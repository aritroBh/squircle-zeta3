import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 1650
t0 = time.time()
from oct5d3_atoms import base_atoms
A = base_atoms()
print("base_atoms at 1650 dps in %.1fs" % (time.time()-t0))
# spot-check against known values
print("Li5(1-b^4) =", mp.nstr(A['Li5(1-b^4)'], 30))
print("Lchi5 =", mp.nstr(A['Lchi5'], 30))
# cross-check Lchi5 vs M4e json (1650 digits)
import json
l5 = json.load(open("<project>/research/08_theory/M4f_w5_extras/Lchi8_5.json"))
print("Lchi5 - M4e 1650-digit:", mp.nstr(abs(A['Lchi5']-mp.mpf(l5['value'])), 3))
# duplication spot-checks: Li5(b^2) ?= 16(Li5(b)+Li5(-b))  -- b^2 not in Li5 list, use Li4(b^2)=16(Li4(b)+Li4(-b))
r2 = mp.sqrt(2); b = r2-1
chk = A['Li4(b^2)'] - 16*(A['Li4(b)']+A['Li4(-b)'])
print("Li4(b^2)-16(Li4(b)+Li4(-b)):", mp.nstr(abs(chk), 3))
chk2 = mp.polylog(5, b**4) - 256*(A['Li5(b)']+A['Li5(-b)'])
print("Li5(b^4)-256(Li5(b)+Li5(-b)):", mp.nstr(abs(chk2), 3))
# Li32 values vs 430-digit OCT5D2 (corrected for their /x bug): theirs*x should equal mine
av = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/atom_vals430.pkl","rb"))
print("Li32(b,1) - b*OCT5D2_L32b1:", mp.nstr(abs(A['Li32(b,1)']-b*mp.mpf(av['L32_b1'])), 3))
print("Li32(h,1) - h*OCT5D2_L32h1:", mp.nstr(abs(A['Li32(1/sqrt2,1)']-(1/r2)*mp.mpf(av['L32_h1'])), 3))
pickle.dump({k: mp.nstr(v, 1630) for k, v in A.items()}, open("atoms_1650.pkl","wb"))
print("saved atoms_1650.pkl")
