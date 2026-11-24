work='<project>/research/08_theory/K2c_kappa2/'
code = r'''
import sympy as sp

sig = sp.Symbol('sig')
q = sp.Symbol('q')
Z = q**2
g1 = sp.Symbol('g1')
s2 = sp.sqrt(2)
ORD = 6

def T(e, n=ORD):
    return sp.series(e, sig, 0, n).removeO()

# saddle for weight e^{-p x}: x_p(x_p+Z) = NZ/p, sig_p = sqrt(p) sig, y_p solves y^2 + q sig_p y = 1
def saddle_p(p):
    sp_ = sp.sqrt(p)*sig
    yp = T(sp.sqrt(1 + q**2*sp_**2/4)) - q*sp_/2
    yp = T(sp.expand(yp))
    xp = T(sp.expand(q*yp/sp_))
    return sp_, yp, xp

# phi_p(x_p) = p x_p + N log(1+Z/x_p) = p q y_p/sig_p - 2 sig^{-2} log y_p
def phi_p(p):
    sp_, yp, xp = saddle_p(p)
    return T(sp.expand(p*q*yp/sp_ - 2*sp.log(1+(yp-1))/sig**2))

# phi_p''(x_p) = N(1/x_p^2 - 1/(x_p+Z)^2)
def phi2_p(p):
    sp_, yp, xp = saddle_p(p)
    return T(sp.expand(sig**-2*((xp+Z)**(-2) - xp**(-2))))*(-1)

D = T(sp.expand(phi_p(2) - phi_p(1)))   # = 2(sqrt2-1)q/sig - Z/2 + c q^3 sig + ...
lead = 2*(s2-1)*q/sig
Dres = T(sp.expand(D - lead))
print("phi2 - phi1 - 2(sqrt2-1)q/sig =", sp.expand(Dres))
E = T(sp.exp(T(sp.expand(-D + lead - sp.Rational(1,2)*Z))))  # e^{-(D-lead)+Z/2}... careful sign
# e^{-(phi2-phi)} = e^{-lead} * e^{-(D-lead)} ; D-lead = -Z/2 + ... so factor = e^{+Z/2}*series
S_exp = T(sp.exp(T(sp.expand(-(D - lead)))))   # = e^{Z/2}(1 + d_exp sig + ...)
S_exp = T(sp.expand(S_exp))
print("exp factor =", sp.expand(S_exp))
# prefactor: [phi2_1/phi2_2]^{1/2} * (x1/x2) * (w2/w1)^{3/2}, w_p = x_p/y_p^2
sp1, y1, x1 = saddle_p(1)
sp2, y2, x2 = saddle_p(2)
w1 = T(sp.expand(x1/y1**2)); w2 = T(sp.expand(x2/y2**2))
P1 = T(sp.expand((phi2_p(1)/phi2_p(2))**(sp.Rational(1,2))))
P2 = T(sp.expand(x1/x2))
P3 = T(sp.expand((w2/w1)**(sp.Rational(3,2))))
Pref = T(sp.expand(P1*P2*P3))
print("prefactor =", sp.expand(Pref))
# bracket ratio: B0(sig_p) = 1 + sig_p/(2q) + O(sig^2)
Br = T(sp.expand((1 + sp2/(2*q))/(1 + sp1/(2*q))))
# m=1/m=0 ratio factor: (1 + (3/2) g1 q sig2/y2')/(1 + (3/2) g1 q sig1/y1') ; m1 weight ~ (3/2) g1 w_p * (1/N) = (3/2) g1 q sig_p/y_p
Mr = T(sp.expand((1 + sp.Rational(3,2)*g1*q*sp2/y2)/(1 + sp.Rational(3,2)*g1*q*sp1/y1)))
R = T(sp.expand((-2)*S_exp*Pref*Br*Mr))
R = T(sp.expand(-R/2))   # divide... wait: total = (-2)*pref(=1/2)*... ; recompute below
# Actually total relative image = (-2)*S_exp*Pref*Br*Mr; S_exp = e^{Z/2}(1+..), Pref = 1/2(1+..)
tot = T(sp.expand(sp.Rational(-2)*S_exp*Pref*Br*Mr))
tot = T(sp.expand(tot/sp.exp(Z/2)))   # divide by e^{Z/2}: expect -(1 + d1 sig + ...)
print("R * e^{+(phi2-phi)_lead} / e^{Z/2} =", sp.expand(tot))
d1 = sp.expand(tot).coeff(sig, 1)
print("d1 =", sp.simplify(d1))
'''
open(work+'step_28.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_28.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-800:])
