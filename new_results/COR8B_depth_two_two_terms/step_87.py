
import mpmath as mp
mp.mp.dps = 40
z = mp.mpf(3)
# Li4(z)+Li4(-z) = (1/8)Li4(z^2)?
lhs = mp.polylog(4,z) + mp.polylog(4,-z)
rhs = mp.polylog(4, z**2)/8
print("Li4(3)+Li4(-3) - Li4(9)/8 =", mp.nstr(lhs-rhs, 12), "(mpmath values)")
lhs2 = mp.conj(mp.polylog(4,z)) + mp.polylog(4,-z)
rhs2 = mp.conj(mp.polylog(4, z**2))/8
print("same with conjugate branch:", mp.nstr(lhs2-rhs2, 12))
# which branch satisfies it?
