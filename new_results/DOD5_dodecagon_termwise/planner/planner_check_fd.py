# Planner check of DOD5's P'''(1) for the regular 12-gon K_1 = {sum_{j=0..5} |<m_j,x>| <= 1}, m_j at angles j*pi/6,
# by polar central finite differences + Richardson (independent of DOD1/DOD2's engines). Also checks the edge length
# L = sqrt(26-15 sqrt3) and the local zeta(3) constant C_6.
import mpmath as mp
mp.mp.dps = 60
ANG = [j*mp.pi/6 for j in range(6)]
N = [(mp.cos(a), mp.sin(a)) for a in ANG]
def P(p):
    def g(th):
        s, c = mp.sin(th), mp.cos(th); F = 0; Fp = 0
        for (a, b) in N:
            v = a*c + b*s; dv = -a*s + b*c
            F += abs(v)**p; Fp += p*abs(v)**(p-1)*mp.sign(v)*dv
        rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*Fp
        return mp.sqrt(rho**2 + drho**2)
    # kinks where some <m_j, (cos th, sin th)> = 0: th = a_j + pi/2 (mod pi); use one period [0, pi/6] by symmetry (12 copies)
    br = [mp.mpf(0), mp.pi/6]   # a_0 + pi/2 = pi/2 = 3*(pi/6): kinks are at multiples of pi/6, so [0, pi/6] is kink-free inside
    return 12*mp.quad(g, br)
def d3(h):
    f = {j: P(1 + j*h) for j in (-2, -1, 1, 2)}
    return (-f[-2] + 2*f[-1] - 2*f[1] + f[2])/(2*h**3)
D = [d3(h) for h in (mp.mpf('2e-4'), mp.mpf('1e-4'), mp.mpf('5e-5'))]
R1 = [(4*D[i+1] - D[i])/3 for i in range(2)]; R2 = (16*R1[1] - R1[0])/15
TGT = mp.mpf('-5.20324042951766621029582970427243461666280533')
print("P(1) =", mp.nstr(P(1), 20), " vs 12*L =", mp.nstr(12*mp.sqrt(26-15*mp.sqrt(3)), 20))
print("raw:", [mp.nstr(x, 12) for x in D]); print("Richardson P'''(1) =", mp.nstr(R2, 15), " DOD5:", mp.nstr(TGT, 15), " diff:", mp.nstr(R2-TGT, 3))
m = 6; s = mp.sin(mp.pi/(2*m)); c = mp.cos(mp.pi/(2*m))
C6 = 48*m*s**4*c*(3*s**2-1); L = mp.sqrt(26-15*mp.sqrt(3))
print("C_6 (paper local law) =", mp.nstr(C6, 25), " (9/4)(2-3sqrt3)L =", mp.nstr(mp.mpf(9)/4*(2-3*mp.sqrt(3))*L, 25))
