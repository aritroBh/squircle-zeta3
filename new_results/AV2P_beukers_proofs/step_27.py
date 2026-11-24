from mpmath import mp, mpf, zeta, log, pi, nprint, catalan, clsin, quad, sqrt
mp.dps = 60
G = catalan; Cl = clsin(2, pi/3)
odd = mpf(open('odd2.txt').read())

W = 2*quad(lambda t: mp.asin(t/2)/(t*(1+t)), [0, 1])
print("W =", nprint(W, 40))
print("W - (8G/3 - pi*log2/3 - Cl) =", nprint(W - (8*G/3 - pi*log(2)/3 - Cl), 6))
print("W - (8G/3 - 2*pi*log2/3 - Cl) =", nprint(W - (8*G/3 - 2*pi*log(2)/3 - Cl), 6))

# odd formula: odd = (pi/8)(2 pi ln2/3 + Cl2 + W)?
candA = pi/8*(2*pi*log(2)/3 + Cl + W)
print("odd - (pi/8)(2 pi ln2/3 + Cl + W):", nprint(odd - candA, 6))
# so required W for odd=pi G/3:
Wreq = 8*G/3 - 2*pi*log(2)/3 - Cl
print("W - Wreq:", nprint(W - Wreq, 6))
# check pieces of W derivation:
I1 = quad(lambda t: log(t)/sqrt(4-t**2), [0,1])
print("int log t/sqrt(4-t^2) =", nprint(I1, 30), " -Cl/2 =", nprint(-Cl/2, 30))
I2 = quad(lambda t: log(1+t)/sqrt(4-t**2), [0,1])
print("int log(1+t)/sqrt(4-t^2) =", nprint(I2, 30), " 4G/3-Cl =", nprint(4*G/3-Cl, 30))
