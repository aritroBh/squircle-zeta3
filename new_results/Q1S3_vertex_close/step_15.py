import mpmath as mp
mp.mp.dps = 50
# check polylog for |z|>1
for z in [mp.mpf('2.4142135623730951'), mp.mpc(1, mp.sqrt(2)), -mp.mpf('2.4142135623730951')]:
    v = mp.polylog(3, z)
    print(z, "->", mp.nstr(v, 20))
# cross-check Li3(-a) via inversion functional equation:
# Li3(z) + Li3(1/z) = -zeta(3) - (1/6)log(-z)^3 - ... standard: Li_n(z)+(-1)^n Li_n(1/z) = -(2 pi i)^n/n! B_n(1/2 + log(-z)/(2 pi i))
# for n=3: Li3(z) - Li3(1/z) = -(2 pi i)^3/6 B_3(1/2 + log(-z)/(2 pi i))
z = -mp.mpf('2.4142135623730951')
w = mp.log(-z)/(2*mp.pi*1j) + mp.mpf('0.5')
B3 = w**3 - w*mp.mpf('0.25')*3/2*2  # B3(x) = x^3 - 3/2 x^2 + 1/2 x
B3 = w**3 - mp.mpf('1.5')*w**2 + mp.mpf('0.5')*w
rhs = -(2*mp.pi*1j)**3/6*B3
print("check:", mp.nstr(mp.polylog(3,z) - mp.polylog(3,1/z) - rhs, 8))
