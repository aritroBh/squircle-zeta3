
import mpmath as mp
mp.mp.dps = 60

def S(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho(phi, m, p):
    return S(phi, m, p)**(-mp.mpf(1)/p)
def speed(phi, m, p):
    r = rho(phi, m, p)
    dr = mp.diff(lambda ph: rho(ph, m, p), phi)
    return mp.sqrt(r*r + dr*dr)
def P_of_p(m, p):
    pts = sorted(set([mp.pi*kk/m for kk in range(2*m)] + [mp.pi/2 + mp.pi*kk/m for kk in range(2*m)]))
    pts = [x for x in pts if 0 < x < 2*mp.pi]
    b = [mp.mpf(0)] + pts + [2*mp.pi]
    return mp.fsum(mp.quad(lambda ph: speed(ph, m, p), [a, bb]) for a, bb in zip(b[:-1], b[1:]))

for m in [2,3]:
    print(f"m={m}: P(1) = {mp.nstr(P_of_p(m, mp.mpf(1)), 30)}  expect {mp.nstr(4*m*(1/mp.cos(mp.pi/(2*m))-mp.cos(mp.pi/(2*m))),30)}")

# finite difference P1 for m=3
for m in [3]:
    for d in [mp.mpf('1e-2'), mp.mpf('1e-3')]:
        fd = (P_of_p(m, 1+d) - P_of_p(m, 1))/d
        print(f"m={m} d={mp.nstr(d,3)}: finite-diff P1 ~ {mp.nstr(fd, 20)}")
    print("paper P1(1) = 2 sqrt3 ln2 =", mp.nstr(2*mp.sqrt(3)*mp.log(2), 20))
    # my derivative integrand approach, but with p-derivatives done analytically:
    # rho(p) = S^(-1/p); log rho = -(1/p) log S
    # d rho/dp = rho * (log S / p^2 - S_p/(p S)), S_p = dS/dp = sum |cos|^p log|cos|
    def dS_dp(phi, m, p):
        out = mp.mpf(0)
        for i in range(1, m+1):
            c = abs(mp.cos(phi - i*mp.pi/m))
            if c > 0:
                out += c**p * mp.log(c)
        return out
    def rho_p1(phi, m):
        Sp = S(phi, m, 1)
        return Sp**(-1) * (mp.log(Sp) - dS_dp(phi, m, 1)/Sp)
    # drho/dphi at p=1, and its p-derivative via diff in phi of rho_p1
    def speed_p1(phi, m):
        r = rho(phi, m, 1)
        dr = mp.diff(lambda ph: rho(ph, m, 1), phi)
        rp = rho_p1(phi, m)
        drp = mp.diff(rho_p1, phi)
        return (r*rp + dr*drp)/mp.sqrt(r*r + dr*dr)
    pts = sorted(set([mp.pi*kk/m for kk in range(2*m)] + [mp.pi/2 + mp.pi*kk/m for kk in range(2*m)]))
    pts = [x for x in pts if 0 < x < 2*mp.pi]
    b = [mp.mpf(0)] + pts + [2*mp.pi]
    P1 = mp.fsum(mp.quad(lambda ph: speed_p1(ph, m), [a, bb]) for a, bb in zip(b[:-1], b[1:]))
    print("analytic-integrand P1(1) =", mp.nstr(P1, 25))
