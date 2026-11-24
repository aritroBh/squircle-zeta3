import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi); a = mp.pi/(2*sqrt2)
b = mp.e**(1j*mp.pi/4)

def Phi1(x, N):   # Phi'(x) = (-1+i)a + N(1/x - b/(1+bx))
    return (-1+1j)*a + N*(1/x - b/(1+b*x))
def Gh(eta, rho, y):
    t = rho*y
    return rho*(eta+t*sqrt2)/(sqrt2*eta*(eta*eta+sqrt2*eta*t+t*t))

random.seed(7)
worst = mp.mpf('0'); warg=None
for _ in range(600):
    y = mp.mpf(10)**random.uniform(-8, -0.4)   # N from ~2 to 1e16
    N = mp.pi/(2*y*y)
    rho = 1 + y*y*random.random()              # true-saddle range; identity is algebraic anyway
    eta = mp.mpf(10)**random.uniform(-1, 1.5)
    x0 = 1/(rho*y)
    lhs = -mp.re(Phi1(eta*x0, N))
    rhs = a - (mp.pi*rho/2)*Gh(eta, rho, y)
    r = abs(lhs-rhs)
    if r > worst: worst, warg = r, (y, rho, eta)
print("L2 identity: worst |lhs-rhs| over 600 random (y,rho,eta):", mp.nstr(worst,6), "at", [mp.nstr(v,4) for v in warg])
# also outside saddle range (arbitrary rho) since the code's identity is pure algebra
worst2 = mp.mpf('0')
for _ in range(300):
    y = mp.mpf(10)**random.uniform(-8, -0.4); N = mp.pi/(2*y*y)
    rho = mp.mpf(10)**random.uniform(-1, 1); eta = mp.mpf(10)**random.uniform(-1, 1.5)
    r = abs(-mp.re(Phi1(eta/(rho*y), N)) - (a - (mp.pi*rho/2)*Gh(eta, rho, y)))
    if r > worst2: worst2 = r
print("L2 identity (arbitrary rho): worst residual:", mp.nstr(worst2,6))

# ---- L3: Gh decreasing in eta ----
import sympy as sp
et, tt, rr = sp.symbols('eta t rho', positive=True)
G = rr*(et+tt*sp.sqrt(2))/(sp.sqrt(2)*et*(et**2+sp.sqrt(2)*et*tt+tt**2))
dG = sp.diff(G, et)
num_target = -(rr/sp.sqrt(2))*(2*et**3 + 4*sp.sqrt(2)*et**2*tt + 4*et*tt**2 + sp.sqrt(2)*tt**3)
den_g = et*(et**2+sp.sqrt(2)*et*tt+tt**2)
diff = sp.simplify(dG - num_target/den_g**2)
print("sympy: Gh'(eta) - [-(rho/rt2)(2e^3+4rt2 e^2 t+4 e t^2+rt2 t^3)/(e(e^2+rt2 e t+t^2))^2] =", diff)
print("=> Gh' numerator is manifestly NEGATIVE for eta,t,rho>0 (all terms positive in -Gh')")
# numeric derivative sign + monotonicity on grid
worst_mono = mp.mpf('0')
for _ in range(2000):
    y = mp.mpf(10)**random.uniform(-8, -0.4); rho = 1+y*y*random.random()
    e1 = mp.mpf(10)**random.uniform(-1, 1.5); e2 = e1*(1+mp.mpf(10)**random.uniform(-6,-0.1))
    v = Gh(e2,rho,y)-Gh(e1,rho,y)   # should be < 0
    if v > worst_mono: worst_mono = v
print("L3 numeric: max Gh(e2)-Gh(e1) over 2000 pairs e2>e1 (should be <0):", mp.nstr(worst_mono,6))
