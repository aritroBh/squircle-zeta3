import mpmath as mp
mp.mp.dps = 50

def S_p(theta, p, m):
    A = mp.pi/m
    return mp.fsum([abs(mp.sin(theta + i*A))**p for i in range(m)])

def Sp_p(theta, p, m):
    A = mp.pi/m
    return mp.fsum([p*abs(mp.sin(theta + i*A))**(p-1)*mp.sign(mp.sin(theta+i*A))*mp.cos(theta + i*A) for i in range(m)])

def rho(theta, p, m):
    return S_p(theta, p, m)**(-mp.mpf(1)/p)

def rho_prime(theta, p, m):
    S = S_p(theta, p, m); Sp = Sp_p(theta, p, m)
    return -S**(-mp.mpf(1)/p - 1)*Sp/p

def density(theta, p, m):
    r = rho(theta, p, m); rp = rho_prime(theta, p, m)
    return mp.sqrt(r*r + rp*rp)

def P_brute2(p, m):
    """Sum of per-segment tanh-sinh integrals over [kA,(k+1)A]."""
    A = mp.pi/m
    tot = mp.mpf(0)
    for k in range(2*m):
        tot += mp.quad(lambda th: density(th, p, m), [k*A, (k+1)*A], method='tanh-sinh')
    return tot

def P_corner_arc(p, m):
    A = mp.pi/m
    return 4*m*mp.quad(lambda th: density(th, p, m), [0, A/2], method='tanh-sinh')

def P_polygonal(p, m, N=20000):
    """Independent: sample boundary densely, sum chord lengths (algebraic accuracy ~1/N^2 per smooth arc)."""
    A = mp.pi/m
    tot = mp.mpf(0)
    prev = None
    for k in range(2*m):
        for j in range(N+1):
            if prev is None and j == 0: pass
            th = k*A + A*j/N
            r = rho(th, p, m)
            pt = (r*mp.cos(th+mp.pi/2), r*mp.sin(th+mp.pi/2))  # arbitrary global frame; chord lengths invariant
            if prev is not None:
                tot += mp.sqrt((pt[0]-prev[0])**2 + (pt[1]-prev[1])**2)
            prev = pt
    return tot

m = 3
for p in [mp.mpf('1.1'), mp.mpf('1.3'), mp.mpf('3')]:
    Pa = P_corner_arc(p, m)
    Pb = P_brute2(p, m)
    Pp = P_polygonal(p, m, 4000)
    print(f"p={mp.nstr(p,5)}: corner-arc={mp.nstr(Pa,30)}")
    print(f"          per-seg  ={mp.nstr(Pb,30)}  res={mp.nstr(abs(Pa-Pb),4)}")
    print(f"          polygonal={mp.nstr(Pp,12)} (low-precision cross-check)")
