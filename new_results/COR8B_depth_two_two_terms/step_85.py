
import mpmath as mp
mp.mp.dps = 40
# branch check
print("Li2(2)  =", mp.nstr(mp.polylog(2,2), 35))
print("expected from above: pi^2/6 + (pi^2/12 - ln2^2/2 - int_1^2 ln(x-1)/x dx) + i*pi*ln2")
# known: Li2(2) = pi^2/3 - (1/2)*ln(2)^2 - i*pi*ln(2) (one convention) or +i*pi*ln2 (other)
z2 = mp.polylog(2,2)
print("pi^2/3 - 0.5*ln2^2 - i*pi*ln2 =", mp.nstr(mp.pi**2/3 - mp.log(2)**2/2 - 1j*mp.pi*mp.log(2), 35))
print("pi^2/3 - 0.5*ln2^2 + i*pi*ln2 =", mp.nstr(mp.pi**2/3 - mp.log(2)**2/2 + 1j*mp.pi*mp.log(2), 35))
# shuffle consistency test: Li2(z)+Li2(1/z) = -pi^2/6 - (1/2) log(-z)^2  for z>1
z = mp.mpf(3)
lhs = mp.polylog(2,z) + mp.polylog(2,1/z)
for sgn in [1,-1]:
    rhs = -mp.pi**2/6 - (mp.log(sgn*-z))**2/2
    print(f"Li2(3)+Li2(1/3) vs -pi^2/6 - log({sgn}*-3)^2/2:", mp.nstr(lhs-rhs, 10))
print("Li2(3)+Li2(1/3) =", mp.nstr(lhs, 35))
