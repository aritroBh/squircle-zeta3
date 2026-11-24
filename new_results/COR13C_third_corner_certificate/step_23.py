import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
t0 = time.monotonic()
from cor13c_common import target2, target3, target4, cross_table, MAPS
from fractions import Fraction as Fr
print("12 maps:", MAPS)
tab = cross_table()
need = ['1/3','-1/3','2/3','1/2','-1/2','-2']
from fractions import Fraction as F
for z in need:
    zz = F(z)
    print(f"cross-ratio {z}: {len(tab.get(zz, []))} reps, first {tab.get(zz, [None])[0]}")
T2 = target2(); print("T2 support:", len(T2), f"({time.monotonic()-t0:.1f}s)")
T3 = target3(); print("T3 support:", len(T3), f"({time.monotonic()-t0:.1f}s)")
T4 = target4(); print("T4 support:", len(T4), f"({time.monotonic()-t0:.1f}s)")
with open('progress.log','a') as f:
    f.write(f"step_23: cor13c_common.py written; target supports {len(T2)}/{len(T3)}/{len(T4)}\n")
