import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.log(2)
H = sp.harmonic
def E_2m(): return H(m)-H(2*m)
def E_2m1(): return H(m-1)-H(2*m-2)+1/(2*m-1)
print("E(2) numeric:", float(E_2m().subs(m,1)), " want:", 1-1/2)
print("E(1) numeric:", float(E_2m1().subs(m,1)), " want:", 1.0)
def Ktval(k, Ek_1):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*Ek_1 - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S
def B_(mm, k): return (-1)**k*pi/(2*mm) - pi*L/2
def Tv(mm, k, Ek_1): return -2*k*Ktval(k, Ek_1) - pi/2 + pi*L/2 - pi/(4*k)
v2 = (pi/2 + B_(4*m+2, 2*m+1) + Tv(4*m+2, 2*m+1, E_2m()))/(4*m+2)
print("v2 at m=1:", sp.simplify(v2.subs(m,1)))
print("want: -13pi/36 + pi ln2/2 =", sp.N(-13*sp.pi/36 + sp.pi*sp.log(2)/2, 15))
print("Ktval(3, E(2)):", sp.simplify(Ktval(3, E_2m().subs(m,1))))
# Kt(6): 2*(K(5)-K(3)+K(1)) = 2*(-pi/48 + pi/16 + pi/8 - pi ln2/4)
print("want Kt(6):", sp.simplify(2*(-sp.pi/48 + sp.pi/16 + sp.pi/8 - sp.pi*sp.log(2)/4)))
