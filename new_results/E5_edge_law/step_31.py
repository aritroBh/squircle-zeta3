import sympy as sp
al, be, u, v, L, h, s0, s, c = sp.symbols('alpha beta u v L h s0 s c', real=True)
w1 = al*u + be*v; e = al*be*u*v; q2m = al**2*u**2 + be**2*v**2
q2 = al**2 + be**2; p = al*be; s1 = al + be; r = h**2/L**2
A5  = 20*w1**3 - 20*e*w1 - 32*q2m - 48*e + 15*w1 - 2
B5  = -15*q2*w1 + 15*p*w1 + 5*p*s1 + 8*q2 - 12*p
A23 = -60*e*w1 + 12*q2m + 72*e - 12*w1 + 2
B23 = -45*p*w1 + 15*p*s1 - 3*q2 + 18*p
Z5law  = 360*al*be*h**2/L*(A5 + r*B5)
Z23law = 120*al*be*h**2/L*(A23 + r*B23)
sub = [(al, 2*s**2), (be, 2*s**2), (u, sp.Rational(1,2)), (v, sp.Rational(1,2)),
       (h, s), (L, 2*s**2/c)]
def redc(expr):  # reduce powers of c: c^2 -> 1-s^2, fully
    expr = sp.expand(expr)
    for pw in range(12, 1, -1):
        expr = sp.expand(expr.subs(c**pw, c**(pw-2)*(1-s**2)))
    return expr
Z5e  = redc(Z5law.subs(sub));  Z23e = redc(Z23law.subs(sub))
tgt5e  = 1440*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)
tgt23e = -480*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)
print("per-edge Z5 - tgt :", sp.factor(redc(Z5e - tgt5e)))
print("per-edge Z23 - tgt:", sp.factor(redc(Z23e - tgt23e)))
