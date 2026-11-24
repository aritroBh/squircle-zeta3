import mpmath as mp
mp.mp.dps = 100

def setup(m):
    A = mp.pi/m
    def pieces(th):
        ls = [mp.log(mp.sin(th + i*A)) for i in range(m)]
        sv = [mp.sin(th + i*A) for i in range(m)]
        cv = [mp.cos(th + i*A) for i in range(m)]
        T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(3)]
        Tp = [mp.fsum(cv)] + [mp.fsum([cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,3)]
        S0 = T[0]; S0p = Tp[0]
        U1 = T[1]/S0
        U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
        U1p = Tp[1]/S0 - T[1]*S0p/S0**2
        lS0 = mp.log(S0); L = S0p/S0
        R1 = lS0 - U1; R2 = -lS0 + U1 - U2
        Q0 = -L; Q1 = L - U1p
        Phi = mp.sqrt(1+Q0**2)
        D0 = Phi/S0
        return R1, R2, Q1, D0, Phi, S0
    return pieces

def J2_integral(m):
    A = mp.pi/m; x = A/2
    pieces = setup(m)
    def integrand(th):
        R1, R2, Q1, D0, Phi, S0 = pieces(th)
        return D0*(R1**2 + 2*R2) + (1/S0)/Phi**3 * Q1**2
    return 4*m*mp.quad(integrand, [0, x], method='tanh-sinh')

def boundary2(m):
    # P(1) * [R1(0)^2 + 2 R2(0)]
    A = mp.pi/m
    P1 = 4*m*mp.sin(A/2)**2/mp.cos(A/2)
    R1, R2, Q1, D0, Phi, S0 = setup(m)(mp.mpf('1e-50'))  # theta->0 limit: R1, R2 finite
    # careful: evaluate limit analytically at 0 instead:
    V = mp.fsum([mp.sin(i*A) for i in range(1,m)])
    T1 = mp.fsum([mp.sin(i*A)*mp.log(mp.sin(i*A)) for i in range(1,m)])
    T2 = mp.fsum([mp.sin(i*A)*mp.log(mp.sin(i*A))**2 for i in range(1,m)])
    S00 = V
    U1 = T1/S00; U2 = T2/(2*S00) - T1**2/(2*S00**2)
    R1 = mp.log(S00) - U1
    R2 = -mp.log(S00) + U1 - U2
    return P1*(R1**2 + 2*R2)

P2_vals = {2: mp.pi**2*mp.sqrt(2)/3,
           3: mp.mpf('-1.00109443973840992962135105887535450171338492998892655629450673110489944312676081385506309'),
           4: mp.mpf('-1.4562051103812046446537306966795936456777152413012283082812637083293482631402138804126453990677')}
for m in [2,3,4]:
    J2 = J2_integral(m)
    bnd = boundary2(m)
    tot = J2 + bnd
    print(f"m={m}: J2 = {mp.nstr(J2,50)}")
    print(f"      boundary = {mp.nstr(bnd,50)},  J2+bnd - P'' = {mp.nstr(abs(tot-P2_vals[m]),5)}")
    P1 = 4*m*mp.sin(mp.pi/(2*m))**2/mp.cos(mp.pi/(2*m))
    print(f"      J2/(P(1)*pi^2) = {mp.nstr(J2/(P1*mp.pi**2),30)}")
