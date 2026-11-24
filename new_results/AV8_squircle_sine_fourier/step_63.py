import sympy as sp
m = sp.symbols('m', integer=True, positive=True)
pi, L = sp.pi, sp.Symbol('L')  # L = ln2

def E(k):  # alternating harmonic sum
    return sp.summation((-1)**(j-1)/j, (j,1,k))

j = sp.symbols('j', integer=True, positive=True)
def Kt(k):
    K1 = pi/8 - pi*L/4
    S = (pi/4)*(-1)**k*(1 - 2*E(k-1) - (-1)**(k-1)/k)
    return 2*(-1)**(k-1)*K1 + S

def A(m_): return pi/2
def B(m_):
    k = m_//2
    return (-1)**k*pi/(2*m_) - pi*L/2
def Tval(m_):
    k = m_//2
    return -2*k*Kt(k) - pi/2 + pi*L/2 - pi/(4*k)

def bn1(n, case):
    # n = 4m+1 (case 1) or 4m+3 (case 3)
    if case == 1:
        n_ = 4*m+1; km, kp = 2*m, 2*m+1
    else:
        n_ = 4*m+3; km, kp = 2*m+1, 2*m+2
    out = 0
    for mm, kk in [(n_-1, km), (n_+1, kp)]:
        out += (A(mm) + B(mm) + Tval(mm))/mm
    return sp.simplify(out/pi)

for case in [1,3]:
    expr = bn1(None, case)
    expr = sp.simplify(expr)
    print("case n=4m+%d:"%case, expr)
