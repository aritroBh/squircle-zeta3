from mpmath import mp, mpf, pi, log, quad, sin, cos, tan, nstr
mp.dps = 40
ln2 = log(2)
def A(m): return quad(lambda u: 1-cos(m*u),[0,pi/2])
def B(m): return quad(lambda u: (1-cos(m*u))*log(cos(u)),[0,pi/2])
def T(m): return quad(lambda u: (1-cos(m*u))*tan(u)**2*log(sin(u)),[0,pi/2])
for m in [2,4]:
    print(m, "A:", nstr(A(m),20), "pi/2:", nstr(pi/2,20))
    print(m, "B:", nstr(B(m),20), "formula:", nstr(-pi/2*ln2 + (-1)**(m//2)*pi/(2*m),20))
    print(m, "T:", nstr(T(m),20))
