# Independent check of the b_n'' law at larger n: direct quadrature of the ODE forcing F2(theta) (closed form),
# b_n'' = (4/pi)/(1-n^2) * int_0^{pi/2} F2 sin(n theta).
import pickle, sympy as sp
from mpmath import mp, mpf, pi, sin, cos, tan, log, clsin, quad, nstr, linspace, euler, zeta
from symsum import t, LN2, GAM, ZETA, Hsym
import symsum
mp.dps = 45
K  = lambda x: (clsin(2, 2*x) + clsin(2, pi-2*x))/2          # sum_{k odd} sin(2kx)/k^2
L0 = lambda x: log(tan(x))
s1  = lambda x: -sin(x)*log(sin(x))/2 - cos(x)*K(x)/2
s1p = lambda x: -cos(x)*log(sin(x))/2 - cos(x)/2 + sin(x)*K(x)/2 + cos(x)*L0(x)/2
F2  = lambda x: -2*L0(x)*s1(x) + 2*tan(x)*s1p(x) - L0(x)**2*sin(x) + 2*L0(x)*sin(x)
law = pickle.load(open('b2_law.pkl', 'rb'))
def lawval(n):
    e = n % 4; tv = (n-e)//4
    sub = {t: tv, LN2: log(2), ZETA[2]: zeta(2)}
    for sym in law[e].free_symbols:
        if str(sym).upper().startswith('H'):
            s_, a = int(str(sym)[1]), int(str(sym).split('_')[1][:-1])
            sub[sym] = sum(mpf(1)/mpf(j)**s_ for j in range(1, a*tv+1))
    return mpf(str(sp.N(law[e].subs(sub), 60)))
for n in (19, 21, 31, 33, 61):
    pts = linspace(0, pi/2, 2*n+1)
    q = quad(lambda x: F2(x)*sin(n*x), pts)*4/pi/(1-n*n)
    print(n, nstr(q, 30), ' law - quad =', nstr(lawval(n)-q, 3), flush=True)
