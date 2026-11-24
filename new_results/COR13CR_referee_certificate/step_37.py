import time
import mpmath as mp
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
import reflib13 as R13
from fractions import Fraction as Fr

tab = R13.cross_values13()
need = [Fr(1,3), Fr(-1,3), Fr(2,3), Fr(1,2), Fr(-1,2), Fr(-2)]
log(f"cross-ratio values of {{0,1/3,1,-1,inf}}: {len(tab)} distinct finite values")
print("  values:", sorted(str(z) for z in tab))
for z in need:
    assert z in tab, z
    log(f"  {z}: {len(tab[z])} frames")
assert Fr(1,9) not in tab
log("  1/9 is NOT a cross-ratio (enters only via the Li4 distribution input): confirmed")

# S3 orbit structure stated in COR13C: orbits {−2,3,−1/2,3/2,1/3,2/3}, {−3,4,−1/3,4/3,1/4,3/4}, {−1,2,1/2}
def anharmonic(z):
    # the anharmonic group: z -> z, 1-z, 1/z, 1/(1-z), z/(z-1), (z-1)/z
    return {z, 1-z, 1/z, 1/(1-z), z/(z-1), (z-1)/z}
orbits = []
seen = set()
for z in sorted(tab, key=str):
    if z in seen: continue
    o = anharmonic(z)
    orbits.append(o)
    seen |= o
log(f"anharmonic-group orbits: {len(orbits)}")
for o in orbits: print("   ", sorted(str(z) for z in o))
expected = [{Fr(-2),Fr(3),Fr(-1,2),Fr(3,2),Fr(1,3),Fr(2,3)},
            {Fr(-3),Fr(4),Fr(-1,3),Fr(4,3),Fr(1,4),Fr(3,4)},
            {Fr(-1),Fr(2),Fr(1,2)}]
assert sorted(map(sorted, orbits)) == sorted(map(sorted, expected))
log("orbits == the three S3-orbits stated in COR13C result.md: confirmed")

# Li4(1/9) = 8(Li4(1/3)+Li4(-1/3)) at 200 dps with direct series (my own)
mp.mp.dps = 210
def li_series(s, z):
    z = mp.mpf(z); tot = mp.mpf(0); n = 1; zn = z
    eps = mp.mpf(10)**(-215)
    while True:
        t = zn/mp.mpf(n)**s
        tot += t
        if n > 100 and abs(t)*abs(z)/(1-abs(z)) < eps: break
        n += 1; zn *= z
    return tot
lhs = li_series(4, mp.mpf(1)/9) - 8*(li_series(4, mp.mpf(1)/3) + li_series(4, -mp.mpf(1)/3))
log(f"|Li4(1/9) - 8(Li4(1/3)+Li4(-1/3))| = {mp.nstr(abs(lhs),3)} at 210 dps (my own series)")
assert abs(lhs) < mp.mpf('1e-200')
with open('progress.log','a') as f:
    f.write("referee step(cross-ratios): all 6 Li args are cross-ratios of {0,1/3,1,-1,inf}, 1/9 absent; orbits match; Li4(1/9) dist 1e-200\n")
print("STEP PASS")
