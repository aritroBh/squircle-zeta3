
import mpmath as mp
mp.mp.dps = 70

def Li21(a, c, N=400):
    """Li_{2,1}(a,c) = sum_{n>m>=1} a^n c^m / (n^2 m), convergent for |a|<=1,|ac|<1"""
    a, c = mp.mpf(a), mp.mpf(c)
    tot = mp.mpf(0)
    for n in range(2, N):
        an = a**n
        s = mp.fsum([an * c**m / (mp.mpf(n)**2 * m) for m in range(1, n)])
        tot += s
    return tot

def Li12(a, c, N=300):
    """Li_{1,2}(a,c) = sum_{n>m>=1} a^n c^m / (n m^2)"""
    a, c = mp.mpf(a), mp.mpf(c)
    tot = mp.mpf(0)
    for n in range(2, N):
        an = a**n
        s = mp.fsum([an * c**m / (n * mp.mpf(m)**2) for m in range(1, n)])
        tot += s
    return tot

def Li11(a, c, N=300):
    a, c = mp.mpf(a), mp.mpf(c)
    tot = mp.mpf(0)
    for n in range(2, N):
        an = a**n
        s = mp.fsum([an * c**m / (n * mp.mpf(m)) for m in range(1, n)])
        tot += s
    return tot

def F(a, b):
    """int_0^1 log(t) log(1-b t)/(1-a t) dt  = Li_{2,1}(a, b/a)/a"""
    a, b = mp.mpf(a), mp.mpf(b)
    return mp.quad(lambda t: mp.log(t)*mp.log(1-b*t)/(1-a*t), [0,1])

k1 = mp.mpf(1)/3; k2 = mp.mpf(1)/9
x1, x2 = -k1, -k2
b1, b2 = k1/(1+k1), k2/(1+k2)   # 1/4, 1/10
ell1 = mp.log(1+k1); ell2 = mp.log(1+k2)
print("beta1 =", mp.nstr(b1,10), " beta2 =", mp.nstr(b2,10))

# the four depth-2 constants
Lam1 = Li21(x1, x2/x1)          # Li21(-1/3, 1/3)
Lam2 = Li21(x2, x1/x2)          # Li21(-1/9, 3)
Lam3 = Li21(b1, b2/b1)          # Li21(1/4, 2/5)
Lam4 = Li21(b2, b1/b2)          # Li21(1/10, 5/2)
print("Lam1 = Li21(-1/3, 1/3)  =", mp.nstr(Lam1, 60))
print("Lam2 = Li21(-1/9, 3)    =", mp.nstr(Lam2, 60))
print("Lam3 = Li21(1/4, 2/5)   =", mp.nstr(Lam3, 60))
print("Lam4 = Li21(1/10, 5/2)  =", mp.nstr(Lam4, 60))

# verify series vs integral F(a,b): Lam1 = F(x1, x2)*x1
print("\n-- series vs integral checks --")
d = abs(F(x1, x2) - Lam1/x1); print("F(x1,x2) vs Lam1/x1:", mp.nstr(d,3))
d = abs(F(x2, x1) - Lam2/x2); print("F(x2,x1) vs Lam2/x2:", mp.nstr(d,3))
# H(a,b) = int Y log(1-bt)/(1-at) ;  D3 = H(x1,x2) = (1/(1+k1))[-ell2 Li2(b1)/b1 + Lam3/b1]
def H(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return mp.quad(lambda t: mp.log(1-t)*mp.log(1-b*t)/(1-a*t), [0,1])
L2b1 = mp.polylog(2, b1); L2b2 = mp.polylog(2, b2)
D3 = ( -ell2*L2b1/b1 + Lam3/b1 ) / (1+k1)
D4 = ( -ell1*L2b2/b2 + Lam4/b2 ) / (1+k2)
print("D3 formula vs quad:", mp.nstr(abs(D3 - H(x1,x2)),3))
print("D4 formula vs quad:", mp.nstr(abs(D4 - H(x2,x1)),3))

# D1..D10 formulas vs quadrature
def Dq(f):
    return mp.quad(f, [0,1])
W1 = lambda t: mp.log(1+k1*t); W2 = lambda t: mp.log(1+k2*t)
X = mp.log; Y = lambda t: mp.log(1-t)
D1 = -Lam1/k1; D2 = -Lam2/k2
print("D1 vs quad:", mp.nstr(abs(D1 - Dq(lambda t: X(t)*W2(t)/(1+k1*t))),3))
print("D2 vs quad:", mp.nstr(abs(D2 - Dq(lambda t: X(t)*W1(t)/(1+k2*t))),3))
rho = k2/(k1-k2); R1 = 1+k1
D5 = (mp.log((k1-k2)/k1)*ell1 + mp.polylog(2,-rho) - mp.polylog(2,-rho*R1))/k1
print("D5 vs quad:", mp.nstr(abs(D5 - Dq(lambda t: W2(t)/(1+k1*t))),3))
D6 = (ell1*ell2 - k1*D5)/k2
print("D6 vs quad:", mp.nstr(abs(D6 - Dq(lambda t: W1(t)/(1+k2*t))),3))
u0 = 1 - k1/k2  # -2
Fu = lambda u: mp.log(u)**2*mp.log(1-u0/u) - 2*mp.log(u)*mp.polylog(2,u0/u) - 2*mp.polylog(3,u0/u) + mp.log(u)**3/3
D7 = (Fu(R1) - Fu(1))/k2
print("D7 vs quad:", mp.nstr(abs(D7 - Dq(lambda t: W1(t)**2/(1+k2*t))),3))
u0b = 1 - k2/k1  # 2/3
R2 = 1+k2
Gu = lambda u: mp.log(u)**2*mp.log(1-u0b/u) - 2*mp.log(u)*mp.polylog(2,u0b/u) - 2*mp.polylog(3,u0b/u) + mp.log(u)**3/3
D8 = (Gu(R2) - Gu(1))/k1
print("D8 vs quad:", mp.nstr(abs(D8 - Dq(lambda t: W2(t)**2/(1+k1*t))),3))
D9  = (ell1**2*ell2 - k2*D7)/(2*k1)
D10 = (ell1*ell2**2 - k1*D8)/(2*k2)
print("D9 vs quad:", mp.nstr(abs(D9  - Dq(lambda t: W1(t)*W2(t)/(1+k1*t))),3))
print("D10 vs quad:", mp.nstr(abs(D10 - Dq(lambda t: W1(t)*W2(t)/(1+k2*t))),3))
print("\nD values:", [mp.nstr(v,20) for v in [D1,D2,D3,D4,D5,D6,D7,D8,D9,D10]])
