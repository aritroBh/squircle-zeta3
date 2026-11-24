import mpmath as mp

mp.mp.dps = 60

def S_p(theta, p, m):
    A = mp.pi/m
    return mp.fsum([abs(mp.sin(theta + i*A))**p for i in range(m)])

def Sp_p(theta, p, m):
    # d/dtheta S_p on a segment where all sines > 0
    A = mp.pi/m
    return mp.fsum([p*mp.sin(theta + i*A)**(p-1)*mp.cos(theta + i*A) for i in range(m)])

def rho(theta, p, m):
    return S_p(theta, p, m)**(-mp.mpf(1)/p)

def rho_prime(theta, p, m):
    S = S_p(theta, p, m); Sp = Sp_p(theta, p, m)
    return -S**(-mp.mpf(1)/p - 1)*Sp/p

def density(theta, p, m):
    r = rho(theta, p, m)
    rp = rho_prime(theta, p, m)
    return mp.sqrt(r*r + rp*rp)

def P_corner_arc(p, m):
    """Perimeter from the exact corner-arc representation: 4m * int_0^{A/2}."""
    A = mp.pi/m
    return 4*m*mp.quad(lambda th: density(th, p, m), [0, A/2], method='tanh-sinh')

def P_brute(p, m):
    """Direct full polar perimeter over [0, 2pi], split at corner directions kA."""
    A = mp.pi/m
    pts = [k*A for k in range(2*m+1)]
    return mp.quad(lambda th: density(th, p, m), pts, method='tanh-sinh')

def P1_formula(m):
    return 4*m*mp.sin(mp.pi/(2*m))**2/mp.cos(mp.pi/(2*m))

print("=== (1a) p=1: corner-arc representation vs closed form vs polygon perimeter ===")
for m in [2,3,4,5,6]:
    A = mp.pi/m
    Pc = P_corner_arc(1, m)
    Pf = P1_formula(m)
    # polygon perimeter directly from corner coordinates
    corners = [mp.matrix([mp.tan(A/2)*mp.sin(k*A - A), mp.tan(A/2)*mp.cos(k*A - A)]) for k in range(2*m+1)]
    Ppoly = mp.fsum([mp.sqrt((corners[k+1][0]-corners[k][0])**2 + (corners[k+1][1]-corners[k][1])**2) for k in range(2*m)])
    print(f"m={m}: corner-arc={mp.nstr(Pc,40)}")
    print(f"      formula  ={mp.nstr(Pf,40)}  polygon={mp.nstr(Ppoly,40)}  res={mp.nstr(abs(Pc-Pf),5)}/{mp.nstr(abs(Pc-Ppoly),5)}")

print("=== (1b) p != 1: corner-arc representation vs brute-force polar perimeter ===")
for m in [3,4]:
    for p in [mp.mpf('0.9'), mp.mpf('1.1'), mp.mpf('1.3'), mp.mpf('2'), mp.mpf('3')]:
        Pa = P_corner_arc(p, m)
        Pb = P_brute(p, m)
        res = abs(Pa-Pb)
        dig = int(-mp.log10(res/abs(Pa))) if res > 0 else 99
        print(f"m={m} p={mp.nstr(p,4)}: P={mp.nstr(Pa,42)}  res={mp.nstr(res,4)}  ({dig} digits)")
