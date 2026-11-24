import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.log(2)
H = sp.harmonic
def Eof(k):
    # E(k) = sum_{j=1}^k (-1)^{j-1}/j, k concrete expression in terms of harmonic
    return sp.Piecewise((H(k//2) - H(k), sp.Even(k)), (H((k-1)//2) - H(k-1) + 1/k, True))
# For symbolic m we know parities: E(2m)=H(m)-H(2m); E(2m-1)=H(m-1)-H(2m-2)+1/(2m-1)... 
def E_2m(): return H(m)-H(2*m)
def E_2m1(): return H(m-1)-H(2*m-2)+1/(2*m-1)   # E(2m-1)
def E_2m2(): return H(m-1)-H(2*m-2)             # E(2m-2)
def E_2m3(): return H(m-2)-H(2*m-4)+1/(2*m-3)   # E(2m-3)
def E_2m4(): return H(m-2)-H(2*m-4)             # E(2m-4)

def Ktval(k, Ek_1):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*Ek_1 - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S

def A_(): return pi/2
def B_(mm, k): return (-1)**k*pi/(2*mm) - pi*L/2
def Tv(mm, k, Ek_1): return -2*k*Ktval(k, Ek_1) - pi/2 + pi*L/2 - pi/(4*k)

# case n=4m+1: k_-=2m (even), k_+=2m+1 (odd); need E(2m-1) and E(2m)
v1 = (A_() + B_(4*m, 2*m) + Tv(4*m, 2*m, E_2m1()))/(4*m)
v2 = (A_() + B_(4*m+2, 2*m+1) + Tv(4*m+2, 2*m+1, E_2m()))/(4*m+2)
expr1 = sp.simplify((v1+v2)/pi)
print("case n=4m+1 simplified:", sp.simplify(expr1))
print("target:", 1/(8*m*(2*m+1)**2))
print("difference:", sp.simplify(expr1 - 1/(8*m*(2*m+1)**2)))
