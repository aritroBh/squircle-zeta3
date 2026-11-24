
import sympy as sp

sig = sp.Symbol('sig', positive=True)
Z   = sp.Symbol('Z')
g1, g2 = sp.symbols('g1 g2')
tau = sp.Symbol('tau', real=True)
sqZ = sp.sqrt(Z)
ORD = 8

# ---- saddle ----
y = 1 - sqZ/2*sig + Z/8*sig**2 - Z**2/128*sig**4 + Z**3/1024*sig**6
xs = sp.expand(sqZ*sig**-1*y)
xs_ser = sp.series(xs, sig, 0, ORD).removeO()

# ---- phi_s ----
logy = sp.series(sp.log(1+(y-1)), sig, 0, ORD).removeO()
phi_s = sp.expand(sqZ*sig**-1*y - 2*sig**-2*logy)
phi_s = sp.series(phi_s, sig, 0, ORD).removeO()
print("phi_s =", phi_s)

# e^{-phi_s} = e^{-2 sqZ/sig} * e^{Z/2} * Efactor, Efactor = series in sig
Efac = sp.expand(-phi_s + 2*sqZ/sig - Z/2)
Efac_ser = sp.series(sp.exp(Efac), sig, 0, ORD).removeO()
print("Efac =", Efac_ser)

# ---- outer integral J_m = int e^{-phi} h_m dx, h_m = (x+Z)^{m+3/2}/x ----
# phase derivatives at x_s: phi(x) = x + N(log(x+Z)-log x), N = 1/sig^2
# phi^{(j)}(x) for j>=2: N*(-1)^j*(j-1)!*[(x+Z)^{-j} - x^{-j}]
x = sp.Symbol('x')
Nsig = sig**-2
def phider(j):
    return Nsig*(-1)**j*sp.factorial(j-1)*((x+Z)**(-j) - x**(-j))
phi2 = sp.expand(phider(2).subs(x, xs_ser))
phi2 = sp.series(phi2, sig, 0, ORD).removeO()
print("phi2(x_s) =", phi2)

# loop expansion: bracket B = <exp(-sum_{j>=3} phi_j tau^j/j!) * (1 + sum_{j>=1} hh_j tau^j/j!)>
# Gaussian with variance 1/phi2 ; keep terms whose order <= sig^2 after moments
# orders: phi_j(x_s) ~ sig^{j-1} ; hh_j ~ sig^j ; <tau^{2r}> ~ sig^{-r}
def series_in_tau(expr, order):
    return sp.series(expr, tau, 0, order).removeO()

mm = sp.Symbol('mm')  # m+3/2
# h(x_s+tau)/h(x_s): log h = mm log(x+Z) - log x
logh = mm*sp.log(x+Z) - sp.log(x)
hh = {}
for j in range(1,7):
    dj = sp.diff(logh, x, j).subs(x, xs_ser)
    dj = sp.series(sp.expand(dj), sig, 0, ORD).removeO()
    hh[j] = dj
phi_ = {}
for j in range(3,9):
    dj = sp.expand(phider(j).subs(x, xs_ser))
    phi_[j] = sp.series(dj, sig, 0, ORD).removeO()

# exponent E(tau) = -sum_{j>=3} phi_j tau^j/j!  ; expand exp(E) in tau, coefficients as sig-series
E = sum(-phi_[j]*tau**j/sp.factorial(j) for j in range(3,9))
expE = sp.series(sp.exp(E), tau, 0, 9).removeO()
Hfac = 1 + sum(hh[j]*tau**j/sp.factorial(j) for j in range(1,7))
integr = sp.expand(expE * Hfac)
integr = sp.series(integr, tau, 0, 9).removeO()
integr = sp.expand(integr)

# Gaussian moments <tau^{2r}> = (2r-1)!! / phi2^r
def moment(r):
    return sp.factorial2(2*r-1)/phi2**r
bracket = 0
for pw in range(0,9,2):
    coeff = sp.expand(integr).coeff(tau, pw)
    r = pw//2
    term = sp.expand(coeff*moment(r))
    bracket = bracket + term
bracket = sp.expand(sp.series(bracket, sig, 0, 4).removeO())
print("bracket (loop) =", bracket)
