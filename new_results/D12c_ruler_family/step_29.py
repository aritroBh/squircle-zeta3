import mpmath as mp
import sympy as sp
mp.mp.dps = 80
eps_s, u_s = sp.symbols('eps u', positive=True)

def Km(qval, m, dps=80):
    mp.mp.dps = dps
    q = sp.Rational(qval)
    e = eps_s
    A = sp.exp(e*(sp.log(u_s)-sp.log(1+u_s))) * (sp.exp(sp.log(1+u_s**(q*(1-e)))/q) - 1)
    Rm = sp.expand(sp.series(A, e, 0, m+1).removeO()).coeff(e, m)
    f = sp.lambdify(u_s, Rm, 'mpmath')
    return (-mp.log(2))**(m+1)/mp.factorial(m+1) + mp.quad(lambda u: f(u)/(u*(1+u)), [0, mp.mpf('0.5'), 1], maxdegree=10)

M0 = mp.log(2) - mp.sqrt(2)*mp.log(1+mp.sqrt(2))
M1 = mp.mpf('2.02095276965635835694690816772523888426901720579')/8
# M2, M3 from squircle note: M3 given; M2 need numeric. Theorem M3: M_3 = 0.022770521491786577456730509742327800250183691729217
M3 = mp.mpf('0.022770521491786577456730509742327800250183691729217')
print("K_0(2) vs M0:", mp.nstr(Km(2,0)-M0,3))
print("K_1(2) vs M1:", mp.nstr(Km(2,1)-M1,3))
print("K_3(2) vs M3:", mp.nstr(Km(2,3)-M3,3))

# Also get M_2 numerically for later: P(p)/8 = 1 + sum M_k p^{-k-1}; fit from P at large p
def P_squircle(p, dps=120):
    mp.mp.dps = dps
    p = mp.mpf(p)
    def ds(th):
        dx = -(2/p)*mp.cos(th)**(2/p-1)*mp.sin(th)
        dy = (2/p)*mp.sin(th)**(2/p-1)*mp.cos(th)
        return mp.sqrt(dx*dx+dy*dy)
    return 4*mp.quad(ds, [0, mp.pi/16, mp.pi/8, mp.pi/4, 3*mp.pi/8, 7*mp.pi/16, mp.pi/2], maxdegree=10)

mp.mp.dps = 100
# Fit M2 using known M0,M1,M3: P/8-1-M0/p-M1/p^2-M3/p^4 = M2/p^3 + M4/p^5...
# use two large p values and known M0,M1,M3 to solve M2 (ignore M4+): need p big enough
for p in [mp.mpf(2000), mp.mpf(4000)]:
    r = P_squircle(p)/8 - 1 - M0/p - M1/p**2 - M3/p**4
    print(f"p={p}: M2 approx {mp.nstr(r*p**3, 30)}")
