import sys
sys.path.insert(0, '.')
import importlib, mpmath as mp
mp.mp.dps = 80
import cor7_lib as L
importlib.reload(L)

BASE = '<project>/research/'
cor5 = mp.mpf(open(BASE+'08_theory/COR5_corner_classicality/M3_2pi3_1100.txt').read().split()[0])
M3, W4, W3, W2 = L.M3_weights(mp.mpf(1)/2)
print("M3(c=1/2) formula =", mp.nstr(M3, 40))
print("COR5 stored       =", mp.nstr(cor5, 40))
print("|diff| =", mp.nstr(abs(M3-cor5), 3))
for nm, mine in [("W4", W4), ("W3", W3), ("W2", W2)]:
    v = mp.mpf(open(BASE+f'08_theory/COR5_corner_classicality/hex_{nm}.txt').read().split()[0])
    print(f"{nm}: diff {mp.nstr(abs(mine-v),3)}")

# also check Psi against G1's stored M3 at another angle (pi/4) -- need M3 stored; G1 values files only have M0..M3? check
print()
for line in open(BASE+'10_frontier/G1_zeta3_any_angle/values_pi_4.txt'):
    print(line.strip())
