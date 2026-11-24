
import sympy as sp

sig = sp.Symbol('sig', positive=True)
Z   = sp.Symbol('Z')
tau = sp.Symbol('tau', real=True)
sqZ = sp.sqrt(Z)
ORD = 8

y = 1 - sqZ/2*sig + Z/8*sig**2 - Z**2/128*sig**4 + Z**3/1024*sig**6
xs_ser = sp.series(sp.expand(sqZ*sig**-1*y), sig, 0, ORD).removeO()
logy = sp.series(sp.log(1+(y-1)), sig, 0, ORD).removeO()
phi_s = sp.series(sp.expand(sqZ*sig**-1*y - 2*sig**-2*logy), sig, 0, ORD).removeO()
Efac = sp.expand(-phi_s + 2*sqZ/sig - Z/2)
Efac_ser = sp.series(sp.exp(Efac), sig, 0, ORD).removeO()

x = sp.Symbol('x'); Nsig = sig**-2
def phider(j):
    return Nsig*(-1)**(j-1)*sp.factorial(j-1)*((x+Z)**(-j) - x**(-j))
phi2 = sp.series(sp.expand(phider(2).subs(x, xs_ser)), sig, 0, ORD).removeO()

mm = sp.Symbol('mm')
# correct Hfac: log h(x_s+tau) - log h(x_s) as series in tau, then exp
lh = mm*(sp.log(x_ser_safe := xs_ser + tau + Z) - sp.log(xs_ser + Z)) - (sp.log(xs_ser + tau) - sp.log(xs_ser))
lh_ser = sp.series(lh, tau, 0, 9).removeO()
lh_ser = sp.series(sp.expand(lh_ser), sig, 0, ORD).removeO()
Hfac = sp.series(sp.exp(lh_ser), tau, 0, 9).removeO()
Hfac = sp.expand(Hfac)

phi_ = {}
for j in range(3,9):
    phi_[j] = sp.series(sp.expand(phider(j).subs(x, xs_ser)), sig, 0, ORD).removeO()
E = sum(-phi_[j]*tau**j/sp.factorial(j) for j in range(3,9))
expE = sp.series(sp.exp(E), tau, 0, 9).removeO()
integr = sp.expand(sp.series(sp.expand(expE*Hfac), tau, 0, 9).removeO())

def moment(r):
    return sp.factorial2(2*r-1)/phi2**r
bracket = 0
for pw in range(0,9,2):
    coeff = integr.coeff(tau, pw)
    bracket += sp.expand(coeff*moment(pw//2))
bracket = sp.expand(sp.series(bracket, sig, 0, 4).removeO())
print("bracket =", bracket)
print()
# prefactor pieces: h_m(x_s) = (NZ)^{m+3/2} x_s^{-(m+5/2)} ; G = sqrt(2pi/phi2)
# h_m(x_s) : use x_s + Z = N Z / x_s = Z sig^{-2}/x_s
hs0 = sp.series(sp.expand((Z*sig**-2)**mm * xs_ser**(-(mm+1))), sig, 0, ORD).removeO()  # mm = m+3/2 -> exponent -(m+5/2) = -(mm+1)
print("h_m(x_s) =", hs0)
G = sp.series(sp.sqrt(2*sp.pi/phi2), sig, 0, ORD).removeO()
print("G =", G)
