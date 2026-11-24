from mpmath import mp, mpf, zeta, log, pi, nprint, sin, quad, pslq, clsin, polylog, catalan
mp.dps = 50
# Honest truncated double sum A_N = sum_{j,k<=N} (1/(jk)) I_jk
def Ijk(j,k):
    if j==k:
        return pi/12 + sin(2*j*pi/3)/(4*j)
    return mpf('0.5')*(sin((j-k)*pi/3)/(2*(j-k)) + sin((j+k)*pi/3)/(2*(j+k)))
A = quad(lambda th: log(2*mp.sin(th))**2, [0, pi/6])
for N in [200, 1000]:
    S = mpf(0)
    for j in range(1, N+1):
        Ijj = pi/12 + sin(2*j*pi/3)/(4*j)
        S += Ijj/j**2
        for k in range(j+1, N+1):
            S += 2*Ijk(j,k)/(j*k)
    print(f"A_N (N={N}) =", nprint(S, 25))
print("A =", nprint(A, 25))
print("A + pi^3/324 =", nprint(A + pi**3/324, 25))
# PSLQ T and Ls3(2pi/3)
T = mpf('1.052682233115594355300917')  # rough; recompute below
def psumH2(N):
    S = mpf(0); H2 = mpf(0)
    for m in range(1, N+1):
        H2 += 1/mpf(m)**2
        S += sin(m*pi/3)*H2/m
    return S
T = psumH2(200000)
Ls3_120 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, 2*pi/3])
G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3)
basis = [1, pi, pi**2, pi**3, z3, L2, L2**2, L2**3, pi**2*L2, pi*L2**2, C13, G,
         polylog(3,mpf('0.5')), pi*C13, C13*L2, pi**2*C13]
names = ['1','pi','pi2','pi3','z3','log2','log2^2','log2^3','pi2log2','pilog2^2','Cl13','G','Li3h','piCl13','Cl13log2','pi2Cl13']
for val, lab in [(T,'T'),(Ls3_120,'Ls3(2pi/3)')]:
    vv = [val]+basis
    r = pslq(vv, tol=mpf('1e-40'), maxcoeff=10**14, maxsteps=10000)
    print(lab, "rel:", r if (r and r[0]!=0) else None)
