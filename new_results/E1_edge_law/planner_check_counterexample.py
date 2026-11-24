# Planner check of the verifier's counterexample: K_p = {|x|^p + |y|^p + 3|x-y|^p <= 4^p}; claimed P'''(1) = 19.9002663631480669804...
import mpmath as mp
mp.mp.dps = 50
def perim(p):
    def f(th):
        c, s = mp.cos(th), mp.sin(th); w = c - s
        S = abs(c)**p + abs(s)**p + 3*abs(w)**p
        dS = p*(abs(c)**(p-1)*mp.sign(c)*(-s) + abs(s)**(p-1)*mp.sign(s)*c + 3*abs(w)**(p-1)*mp.sign(w)*(-s - c))
        rho = 4*S**(-1/p); drho = 4*(-1/p)*S**(-1/p-1)*dS
        return mp.sqrt(rho**2 + drho**2)
    return 2*mp.quad(f, [0, mp.pi/4, mp.pi/2, mp.pi])
print("P(1) own", mp.nstr(perim(1), 25), " exact 4sqrt5+2sqrt2 =", mp.nstr(4*mp.sqrt(5) + 2*mp.sqrt(2), 25))
h = mp.mpf('0.004'); K = 7; xs = [k*h for k in range(-K, K+1)]
ys = [perim(1 + x) for x in xs]
co = mp.lu_solve(mp.matrix([[x**j for j in range(2*K+1)] for x in xs]), mp.matrix(ys))
P3 = 6*co[3]
z3, pi, L = mp.zeta(3), mp.pi, mp.log(2)
closed = mp.sqrt(5)*(mp.mpf(231)/125*z3 + mp.mpf(111)/125*pi**2*L + mp.mpf(52)/125*L**3 - mp.mpf(16)/25*pi**2 - mp.mpf(288)/25*L**2 + mp.mpf(96)/5*L) \
       + mp.sqrt(2)*(-mp.mpf(15)/64*z3 + pi**2*L/16 + 16*L**3 - pi**2/16 - 48*L**2 + 24*L)
print("P'(1) own", mp.nstr(co[1], 15), "  P''(1) own", mp.nstr(2*co[2], 15))
print("P'''(1) own", mp.nstr(P3, 16), "  claimed", mp.nstr(closed, 16), "  diff", mp.nstr(P3 - closed, 3))
