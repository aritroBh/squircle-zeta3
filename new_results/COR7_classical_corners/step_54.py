import sys, json
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L

mp.mp.dps = 1120
c = mp.mpf(1)/3
M3 = mp.mpf(open("c13_M3_hpl_1150.txt").read())
W4v = mp.mpf(open("c13_W4_hpl_1150.txt").read())
W3v = mp.mpf(open("c13_W3_hpl_1150.txt").read())
W2v = mp.mpf(open("c13_W2_hpl_1150.txt").read())

l2, l3 = mp.log(2), mp.log(3)
z2, z3 = mp.zeta(2), mp.zeta(3)
pi2 = mp.pi**2
def Li(s, x): return mp.polylog(s, x)
l313, l31m, l213 = Li(3, c), Li(3, -mp.mpf(1)/2), Li(2, c)

# W2 exact from formula
W2_cf = -Li(2, c) - Li(2, -c) + pi2/12
print("W2 - (-Li2(1/3)-Li2(-1/3)+pi^2/12) =", mp.nstr(abs(W2v - W2_cf), 3))

# W3 closed form
W3_cf = (-26*z3 - 87*l313 - 130*l31m + mp.mpf(21)/4*pi2*l2 - mp.mpf(43)/12*pi2*l3
         + 67*l213*l2 - 45*l213*l3 + mp.mpf(65)/3*l2**3 - 65*l2**2*l3 + mp.mpf(109)/2*l2*l3**2 - mp.mpf(37)/3*l3**3)
print("W3 - closed form =", mp.nstr(abs(W3v - W3_cf), 3))

# W4 closed form (from c13_F4_relation.json)
rel = json.load(open("c13_F4_relation.json"))
a0, coeffs = rel["a0"], rel["coeffs"]
B = [Li(4, -2), Li(4, c), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
     l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
     l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
     pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
     l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
W4_cf = -sum(mp.mpf(ci)*b for ci, b in zip(coeffs, B))/a0
print("W4 - closed form =", mp.nstr(abs(W4v - W4_cf), 3))

# total
print("M3 - (W4cf+W3cf+W2cf) =", mp.nstr(abs(M3 - (W4_cf + W3_cf + W2_cf)), 3))

# Psi^(1)(1/3) closed form: (21/2)log2 - (15/2)log3
P1 = L._combo(L.PSI, c, wts={1})
print("Psi1(1/3) - [(21/2)log2-(15/2)log3] =", mp.nstr(abs(P1 - (mp.mpf(21)/2*l2 - mp.mpf(15)/2*l3)), 3))
P2 = L._combo(L.PSI, c, wts={2})
print("Psi2(1/3) =", mp.nstr(P2, 30))
# report zeta(3) coefficients
print("\npure zeta(3) coefficient in M3 closed form (basis {Li3(1/3), Li3(-1/2)}): -26")
print("zeta(3)*log terms in W4: (14/3) z3 log2 - (25/9) z3 log3")
