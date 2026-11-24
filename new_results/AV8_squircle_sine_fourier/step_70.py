import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.log(2)
H = sp.harmonic
def E_2m(): return H(2*m)-H(m)
def E_2m1(): return H(2*m)-H(m)+1/(2*m)     # E(2m-1)
def E_2mp1(): return H(2*m+2)-H(m+1)         # E(2m+2) for case3
def Ktval(k, Ek_1):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*Ek_1 - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S
def B_(mm, k): return (-1)**k*pi/(2*mm) - pi*L/2
def Tv(mm, k, Ek_1): return -2*k*Ktval(k, Ek_1) - pi/2 + pi*L/2 - pi/(4*k)

# case n=4m+1
v1 = (pi/2 + B_(4*m, 2*m) + Tv(4*m, 2*m, E_2m1()))/(4*m)
v2 = (pi/2 + B_(4*m+2, 2*m+1) + Tv(4*m+2, 2*m+1, E_2m()))/(4*m+2)
expr1 = sp.simplify((v1+v2)/pi)
d1 = sp.simplify(expr1 - 1/(8*m*(2*m+1)**2))
d1 = d1.subs({H(m-1): H(m)-1/m, H(2*m-2): H(2*m)-1/(2*m-1)-1/(2*m), H(2*m-1): H(2*m)-1/(2*m)})
print("case n=4m+1: expr - law =", sp.simplify(d1))
