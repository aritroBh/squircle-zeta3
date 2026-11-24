import sys
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/10_frontier/G1_zeta3_any_angle')
import mpmath as mp
mp.mp.dps = 60
import corner_series as cs

c = mp.mpf(1)/3
alpha = 2*mp.acos(c)   # full interior angle

# M0 closed form from G1: M0(c) = -(1+c)log(1+c) - (1-c)log(1-c), c = cos(alpha/2)
M0 = -(1+c)*mp.log(1+c) - (1-c)*mp.log(1-c)
print("M0 formula (c=cos(a/2)=1/3):", mp.nstr(M0, 40))

for a_try, name in [(alpha, "a=alpha(full)"), (alpha/2, "a=alpha/2")]:
    M0q = cs.corner_coeffs(a_try, 0)[0]
    print(f"quadrature with {name}: {mp.nstr(M0q, 40)}  diff {mp.nstr(abs(M0q-M0),3)}")
