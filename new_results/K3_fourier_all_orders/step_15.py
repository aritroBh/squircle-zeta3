import sympy as sp
import mpmath as mp
import pickle

mp.mp.dps = 60
pi = mp.pi; z = mp.zeta
u = sp.symbols('u')  # u = k^2

# Q_n as polynomial in u=k^2 with zeta-product coefficients, from PSLQ results (step_14)
# Q_n(kappa), kappa^2 = 4u;  q_{n,r} (coeff of kappa^{2r}) = (-1)^r D~_{n,r}/(2r)!
Dt = {  # D~_{n,r} coefficient dicts: zeta monomial -> rational
 (2,0): {('z',2): sp.Rational(1)},
 (3,0): {('z',3): sp.Rational(1)},
 (4,0): {('z',4): sp.Rational(-33,16)},
 (4,1): {('z',4): sp.Rational(7,2)},
 (5,0): {('z',5): sp.Rational(-3,4), ('p',(2,3)): sp.Rational(-7,4)},
 (5,1): {('z',5): sp.Rational(2), ('p',(2,3)): sp.Rational(2)},
 (6,0): {('z',6): sp.Rational(4039,256), ('p',(3,3)): sp.Rational(-7,8)},
 (6,1): {('z',6): sp.Rational(-733,16), ('p',(3,3)): sp.Rational(1)},
 (6,2): {('z',6): sp.Rational(93,2)},
 (7,0): {('z',7): sp.Rational(129,16), ('p',(2,5)): sp.Rational(141,16), ('p',(2,2,3)): sp.Rational(1071,160)},
 (7,1): {('z',7): sp.Rational(-23), ('p',(2,5)): sp.Rational(-25), ('p',(2,2,3)): sp.Rational(-181,10)},
 (7,2): {('z',7): sp.Rational(24), ('p',(2,5)): sp.Rational(24), ('p',(2,2,3)): sp.Rational(84,5)},
 (8,0): {('z',8): sp.Rational(-2891923,12288), ('p',(3,5)): sp.Rational(141,16), ('p',(2,3,3)): sp.Rational(169,32)},
 (8,1): {('z',8): sp.Rational(495965,512), ('p',(3,5)): sp.Rational(-25), ('p',(2,3,3)): sp.Rational(-27,2)},
 (8,2): {('z',8): sp.Rational(-60889,32), ('p',(3,5)): sp.Rational(24), ('p',(2,3,3)): sp.Rational(12)},
 (8,3): {('z',8): sp.Rational(5715,4)},
 (9,0): {('z',9): sp.Rational(-22805,192), ('p',(2,7)): sp.Rational(-8343,64), ('p',(2,2,5)): sp.Rational(-60093,640), ('p',(2,2,2,3)): sp.Rational(-263713,4480), ('p',(3,3,3)): sp.Rational(169,96)},
 (9,1): {('z',9): sp.Rational(3907,8), ('p',(2,7)): sp.Rational(4163,8), ('p',(2,2,5)): sp.Rational(29741,80), ('p',(2,2,2,3)): sp.Rational(130973,560), ('p',(3,3,3)): sp.Rational(-9,2)},
 (9,2): {('z',9): sp.Rational(-958), ('p',(2,7)): sp.Rational(-990), ('p',(2,2,5)): sp.Rational(-3501,5), ('p',(2,2,2,3)): sp.Rational(-3093,7), ('p',(3,3,3)): sp.Rational(4)},
 (9,3): {('z',9): sp.Rational(720), ('p',(2,7)): sp.Rational(720), ('p',(2,2,5)): sp.Rational(504), ('p',(2,2,2,3)): sp.Rational(2232,7)},
}
def monovalue(mo):
    if mo[0]=='z': return z(mo[1])
    return mp.fprod([z(i) for i in mo[1]])
def Q(n, uval):
    tot = mp.mpf(0)
    for r in range((n-2)//2+1):
        q = sp.Rational((-1)**r, sp.factorial(2*r))
        d = Dt[(n,r)]
        tot += mp.mpf(q) * (4*uval)**r * mp.fsum(mp.mpf(str(c))*monovalue(mo) for mo,c in d.items())
    return tot
def cn(n, k):  # c_n(k)
    return (-1)**k * mp.sqrt(2)/pi * Q(n, k*k)

# verify against f1_03 stored c_7 values
print("c_7 check vs f1_03 stored values:")
for k, val in [(0, '20.22684832'), (1, '-168.4470767'), (2, '1044.464587'), (5, '-25293.57381')]:
    print(f"  k={k}: mine {mp.nstr(cn(7,k),12)}  f1 {val}  diff {mp.nstr(abs(cn(7,k)-mp.mpf(val)),3)}")
# verify c_4, c_5 against F1 formulas
print("\nc_4: mine vs -(33/16+7k^2)z4:", [mp.nstr(cn(4,k)/((-1)**k*mp.sqrt(2)/pi) + (sp.Rational(33,16)+7*k*k)*z(4), 3) for k in range(4)])
print("c_5: mine vs -((3+16u)z5+(7+16u)z2z3)/4:", [mp.nstr(cn(5,k)/((-1)**k*mp.sqrt(2)/pi) + ((3+16*k*k)*z(5)+(7+16*k*k)*z(2)*z(3))/4, 3) for k in range(4)])

# pretty print Q_n(u)
zeta = sp.symbols('z2:12')
def Qsym(n):
    res = sp.Rational(0)
    for r in range((n-2)//2+1):
        q = sp.Rational((-1)**r, sp.factorial(2*r))
        for mo, c in Dt[(n,r)].items():
            if mo[0]=='z': sym = zeta[mo[1]-2]
            else: sym = sp.prod([zeta[i-2] for i in mo[1]])
            res += q*c*sym*(4*u)**r
    return sp.expand(res)
print("\nQ_n(u), u=k^2:  c_n(k) = (-1)^k (sqrt2/pi) Q_n(k^2)")
for n in range(2, 10):
    print(f"Q_{n} =", Qsym(n))
with open("Qpolys.pkl","wb") as f:
    pickle.dump({n: str(Qsym(n)) for n in range(2,10)}, f)
