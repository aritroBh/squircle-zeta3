import sys, json
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
import mpmath as mp
from cor7_pt import QuadField, QuadPoint
from cor6_sym import W4Quotient, lneg
from cor6_pipe import W4TERMS
from cor6_ft import ft_coords

Q2 = W4Quotient(2)
def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})

QF = QuadField(2)
P = QuadPoint(QF, "b=sqrt2-1", (Fr(-1), Fr(1)), primes=[2])
F = P.F
rt2 = (Fr(0), Fr(1))
# Li4 args as NF elements
B = lambda: (Fr(-1), Fr(1))
args = {"b": (Fr(-1), Fr(1)), "-b": (Fr(1), Fr(-1)), "-b^2": (Fr(-3), Fr(2)),
        "1/sqrt2": (Fr(0), Fr(1,2)), "1/2": (Fr(1,2), Fr(0)),
        "1-1/sqrt2": (Fr(1), Fr(-1,2)), "2-sqrt2": (Fr(2), Fr(-1))}
# coefficients (p + q sqrt2) from the W4 relation
coef = {"b": (Fr(22), Fr(-31,2)), "-b": (Fr(24), Fr(-27,2)), "-b^2": (Fr(-5,2), Fr(13,8)),
        "1/sqrt2": (Fr(0), Fr(-1)), "1/2": (Fr(3,4), Fr(0)),
        "1-1/sqrt2": (Fr(-4), Fr(4)), "2-sqrt2": (Fr(-2), Fr(3))}
comp0 = [Fr(0)]*3; comp1 = [Fr(0)]*3
for nm, x in args.items():
    v = li4_coords(P, x)
    p, q = coef[nm]
    for j in range(3):
        comp0[j] += p*v[j]
        comp1[j] += q*v[j]
print("Li4 part of numerical W4: theta^0:", tuple(comp0), " theta^1:", tuple(comp1))

# exact symbol class of Omega(b)
SOm = P.omega_symbol(W4TERMS)
cls = ft_coords(Q2, SOm, 2)
print("class(Omega(sqrt2-1)) exact:    theta^0:", cls[0], " theta^1:", cls[1])
print("MATCH:", tuple(comp0) == tuple(cls[0]) and tuple(comp1) == tuple(cls[1]))

# total assembly at dps 1120
mp.mp.dps = 1120
M3 = mp.mpf(open("b21_M3_hpl_1150.txt").read())
b = mp.sqrt(2)-1; rt2v = mp.sqrt(2)
Llog = mp.log(1+rt2v); lam = mp.log(2)
def Li(s, x): return mp.polylog(s, x)
z3, pi2 = mp.zeta(3), mp.pi**2
rel4 = json.load(open("b21_W4_relation.json"))
a0, cf = rel4["a0"], rel4["coeffs"]
Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2v)
Li2b = Li(2, b)
B4 = [Li(4, b), Li(4, -b), Li(4, -b*b), Li(4, 1/rt2v), Li(4, mp.mpf(1)/2), Li(4, 1-1/rt2v), Li(4, 2-rt2v),
      Llog*Li3b, Llog*Li3mb, Llog*Li3s, Llog*z3, lam*Li3b, lam*Li3mb, lam*Li3s, lam*z3,
      Li2b**2, Li2b*Llog**2, Li2b*lam*Llog, Li2b*lam**2, pi2*Li2b,
      mp.pi**4, pi2*Llog**2, pi2*lam*Llog, pi2*lam**2,
      Llog**4, Llog**3*lam, Llog**2*lam**2, Llog*lam**3, lam**4]
nb = len(B4)
W4cf = -sum(mp.mpf(cf[i])*B4[i] + mp.mpf(cf[i+nb])*rt2v*B4[i] for i in range(nb))/a0
rel3 = json.load(open("b21_W3_relation.json"))
a03, cf3 = rel3["a0"], rel3["coeffs"]
B3 = [z3, Li3b, Li3mb, Li3s, Llog*Li2b, lam*Li2b, pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
nb3 = len(B3)
W3cf = -sum(mp.mpf(cf3[i])*B3[i] + mp.mpf(cf3[i+nb3])*rt2v*B3[i] for i in range(nb3))/a03
W2cf = -Li(2, b*b)/2 + pi2/12
print("M3 - (W4cf+W3cf+W2cf) =", mp.nstr(abs(M3 - (W4cf+W3cf+W2cf)), 3))
W4v = mp.mpf(open("b21_W4_hpl_1150.txt").read()); W3v = mp.mpf(open("b21_W3_hpl_1150.txt").read()); W2v = mp.mpf(open("b21_W2_hpl_1150.txt").read())
print("|W4-W4cf| =", mp.nstr(abs(W4v-W4cf),3), " |W3-W3cf| =", mp.nstr(abs(W3v-W3cf),3), " |W2-W2cf| =", mp.nstr(abs(W2v-W2cf),3))
