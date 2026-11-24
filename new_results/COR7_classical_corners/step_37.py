import sys
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/10_frontier/G1_zeta3_any_angle')
import mpmath as mp
mp.mp.dps = 70
import cor7_lib as L
import corner_series as cs

c = mp.mpf(1)/3
alpha = 2*mp.acos(c)
M3, W4, W3, W2 = L.M3_weights(c)
M3q = cs.corner_coeffs(alpha, 3)[3]
print("M3 HPL formula =", mp.nstr(M3, 50))
print("M3 quadrature  =", mp.nstr(M3q, 50))
print("|diff| =", mp.nstr(abs(M3-M3q), 3))
print("W4+W3+W2-M3 =", mp.nstr(abs(W4+W3+W2-M3),3))
print("W4 =", mp.nstr(W4, 30)); print("W3 =", mp.nstr(W3, 30)); print("W2 =", mp.nstr(W2, 30))
