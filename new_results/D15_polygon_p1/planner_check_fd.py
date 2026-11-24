# Planner check of D15: own polar perimeter of K_p = {sum_i |<n_i,x>|^p <= 1}, n_i = e(i pi/m); Taylor fit around p = 1.
import mpmath as mp
mp.mp.dps = 50
def perim(m, p):
    ang = [i*mp.pi/m for i in range(1, m+1)]
    def f(th):
        c = [mp.cos(th - a) for a in ang]; s = [mp.sin(th - a) for a in ang]
        S = mp.fsum(abs(x)**p for x in c)
        Sp = mp.fsum(-p*abs(x)**(p-1)*mp.sign(x)*y for x, y in zip(c, s))
        rho = S**(-1/p); drho = -(1/p)*S**(-1/p-1)*Sp
        return mp.sqrt(rho**2 + drho**2)
    kinks = sorted(set([mp.mpf(0), mp.pi] + [ (a + mp.pi/2) % mp.pi for a in ang]))
    return 2*mp.quad(f, kinks)
for m in (2, 3):
    print("m", m, "P(1)", mp.nstr(perim(m, 1), 20), "claimed", mp.nstr(4*m*(mp.sec(mp.pi/(2*m)) - mp.cos(mp.pi/(2*m))), 20))
m = 3; h = mp.mpf('0.004'); K = 7
xs = [k*h for k in range(-K, K+1)]
ys = [perim(m, 1 + x) for x in xs]
# least-squares-free: exact polynomial interpolation of degree 2K, read off Taylor coefficients
coeffs = mp.lu_solve(mp.matrix([[x**j for j in range(2*K+1)] for x in xs]), mp.matrix(ys))
d1, d2, d3 = coeffs[1], 2*coeffs[2], 6*coeffs[3]
r3, L = mp.sqrt(3), mp.log(2)
cl1 = 2*r3*L
cl2 = r3/8*mp.pi**2 + 2*r3*L**2 - 4*r3*L
cl3 = r3/8*(16*L**3 + 3*mp.pi**2*L - 96*L**2 + 96*L - 9*mp.zeta(3) - 6*mp.pi**2)
for nm, a, b in (("P'", d1, cl1), ("P''", d2, cl2), ("P'''", d3, cl3)):
    print(nm, "fit", mp.nstr(a, 16), "claimed", mp.nstr(b, 16), "rel diff", mp.nstr(abs(a-b)/abs(b), 3))
print("P(1.1) own", mp.nstr(perim(3, mp.mpf('1.1')), 30), "D15 table 3.69840169493827649453206212737602")
