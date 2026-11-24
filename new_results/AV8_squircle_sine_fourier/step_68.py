import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.log(2)
H = sp.harmonic
def E_2m(): return H(m)-H(2*m)
def E_2m1(): return H(m-1)-H(2*m-2)+1/(2*m-1)
def Ktval(k, Ek_1):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*Ek_1 - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S
def B_(mm, k): return (-1)**k*pi/(2*mm) - pi*L/2
def Tv(mm, k, Ek_1): return -2*k*Ktval(k, Ek_1) - pi/2 + pi*L/2 - pi/(4*k)
# check pieces at m=1 (n=5): want (A+B+T)(4) = 3pi/8*4? no: (3/2 pi - 2 pi ln2); /4
v1 = (pi/2 + B_(4*m, 2*m) + Tv(4*m, 2*m, E_2m1()))/(4*m)
print("v1 at m=1:", sp.simplify(v1.subs(m,1)), "  want: 3pi/32 - pi ln2/8? no, (3/2 pi -2 pi L)/4 = 3pi/8 - pi L/2")
print("Ktval(2, E(1)) at m=1:", sp.simplify(Ktval(2, E_2m1().subs(m,1))))
# Kt(2) from Fraction code: 2*K(1) = pi/4 - pi ln2/2
print("want Kt(2): pi/4 - pi ln2/2 =", sp.nsimplify(sp.pi/4 - sp.pi*sp.log(2)/2))
