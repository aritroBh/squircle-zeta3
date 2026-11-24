# Planner check of the general edge law on a NON-rhombus parallelogram K_p = {|x|^p + |x+2y|^p <= 1}:
# predicted P = P(1) + (pi^2 S/6) d^2 + S (zeta3 - pi^2/3) d^3 + O(d^4), S = sum_e h_e^2/L_e = 1/(4 sqrt2) + 1/2.
import mpmath as mp
mp.mp.dps = 50
forms = [(1, 0), (1, 2)]
def perim(p):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        vals = [a*c + b*s for a, b in forms]; ders = [-a*s + b*c for a, b in forms]
        S = mp.fsum(abs(v)**p for v in vals)
        Sp = mp.fsum(p*abs(v)**(p-1)*mp.sign(v)*dv for v, dv in zip(vals, ders))
        rho = S**(-1/p); drho = -(1/p)*S**(-1/p-1)*Sp
        return mp.sqrt(rho**2 + drho**2)
    kinks = sorted({mp.mpf(0), mp.pi, mp.pi/2, mp.atan2(1, -2) % mp.pi})
    return 2*mp.quad(f, kinks)
S = 1/(4*mp.sqrt(2)) + mp.mpf(1)/2
P1 = 2*(mp.sqrt(2) + 1)
print("P(1) own", mp.nstr(perim(1), 20), " vs 2(sqrt2+1)", mp.nstr(P1, 20))
h = mp.mpf('0.004'); K = 7; xs = [k*h for k in range(-K, K+1)]
ys = [perim(1 + x) for x in xs]
co = mp.lu_solve(mp.matrix([[x**j for j in range(2*K+1)] for x in xs]), mp.matrix(ys))
P1d, P2d, P3d = co[1], 2*co[2], 6*co[3]
print("P'(1)  own", mp.nstr(P1d, 12), " predicted 0")
print("P''(1) own", mp.nstr(P2d, 16), " predicted pi^2 S/3 =", mp.nstr(mp.pi**2*S/3, 16))
print("P'''(1) own", mp.nstr(P3d, 14), " predicted 6 S (zeta3 - pi^2/3) =", mp.nstr(6*S*(mp.zeta(3) - mp.pi**2/3), 14))
