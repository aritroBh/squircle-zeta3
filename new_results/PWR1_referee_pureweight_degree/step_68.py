import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 105
def eta(s): return (1-mp.power(2,1-s))*mp.zeta(s)
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m*mp.log1p(mp.e**(-s))**j, [0, 0.5, 2, 8, 32, mp.inf])
    return Icache[(m,j)]
e = {}
with open('e_rm_mine.pkl','rb') as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(',')); e[(r,m)] = mp.mpf(val)
def Dtilde(n, r):
    tot = mp.mpf(0)
    for j in range(1, n-2*r):
        m = n-1-j
        tot += e[(r,m)]*I(m,j)/mp.factorial(j)
    return 2*tot

# Theta_n(u) = sum_m D~_{n,m} (-4u)^m/(2m)!   (my independent build)
def Theta(n, u):
    return mp.fsum(Dtilde(n, m)*(-4*u)**m/mp.factorial(2*m) for m in range(0, n//2))

# K3's stored Qpolys
Q = pickle.load(open('../K3_fourier_all_orders/Qpolys.pkl','rb'))
print(type(Q), str(Q)[:400])
