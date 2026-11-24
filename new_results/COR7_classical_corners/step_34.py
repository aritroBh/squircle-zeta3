import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 80
import cor7_lib as L

BASE = '<project>/research/'
# stored M3 at 2pi/3 from G1 and COR5
stored = None
for line in open(BASE+'10_frontier/G1_zeta3_any_angle/values_2pi_3.txt'):
    if line.startswith('M3 ='):
        stored = mp.mpf(line.split('=')[1].strip())
cor5 = mp.mpf(open(BASE+'08_theory/COR5_corner_classicality/M3_2pi3_1100.txt').read().split()[0])
print("G1 stored M3(2pi/3)  =", mp.nstr(stored, 60))
print("COR5 M3(2pi/3)       =", mp.nstr(cor5, 60))
M3, W4, W3, W2 = L.M3_weights(mp.mpf(1)/2)
print("my M3 formula (c=1/2)=", mp.nstr(M3, 60))
print("|diff G1| =", mp.nstr(abs(M3-stored),3))
# also compare weights with COR5 hex files
for nm, mine in [("W4", W4), ("W3", W3), ("W2", W2)]:
    v = mp.mpf(open(BASE+f'08_theory/COR5_corner_classicality/hex_{nm}.txt').read().split()[0])
    print(f"{nm}: mine {mp.nstr(mine, 25)} vs COR5 {mp.nstr(v, 25)}  diff {mp.nstr(abs(mine-v),3)}")
