import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 105
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
def Theta(n, u):
    return mp.fsum(Dtilde(n, m)*(-4*u)**m/mp.factorial(2*m) for m in range(0, n//2))

Q = pickle.load(open('../K3_fourier_all_orders/Qpolys.pkl','rb'))
u_sym = sp.Symbol('u')
zsubs = {sp.Symbol(f'z{k}'): sp.zeta(k) for k in range(2, 13)}
print("Theta_n(u) [my D~-build, quadrature I's] vs K3's Q_n(u) [PSLQ zeta-products]:")
for n in range(2, 9):
    qn = sp.sympify(Q[n], locals={f'z{k}': sp.Symbol(f'z{k}') for k in range(2, 13)})
    diffs = []
    for uv in [0, 1, 2, 5, 7]:
        a = Theta(n, uv)
        b = mp.mpf(str(qn.subs({u_sym: uv, **zsubs}).evalf(50)))
        diffs.append(abs(a-b))
    print(f"  n={n}: max|diff| over u in {{0,1,2,5,7}} = {mp.nstr(max(diffs), 3)}")
print("\npaper claim: leading coeff of Theta_{2m} = (-1)^{m-1}(2^{2m-1}-1) zeta(2m):")
for mv in range(1, 7):
    n = 2*mv
    top = Dtilde(n, mv-1)*(-4)**(mv-1)/mp.factorial(2*mv-2)
    claim = (-1)**(mv-1)*(2**(2*mv-1)-1)*mp.zeta(2*mv)
    print(f"  m={mv}: top={mp.nstr(top,30)} diff={mp.nstr(abs(top-claim),3)}")
