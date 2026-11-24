# step_05.py: THE ZETA(5) EDGE LAW -- closed forms and symbolic verification.
# Z_e^(5)  = 360 al be h^2/L [ A5 + (h^2/L^2) B5 ],  Z_e^(2,3) = 120 al be h^2/L [ A23 + (h^2/L^2) B23 ]
# with w1=al u+be v, e=al be u v, q2m=al^2 u^2+be^2 v^2, q2=al^2+be^2, p=al be, s1=al+be:
#   A5  = 20 w1^3 - 20 e w1 - 32 q2m - 48 e + 15 w1 - 2
#   B5  = -15 q2 w1 + 15 p w1 + 5 p s1 + 8 q2 - 12 p
#   A23 = -60 e w1 + 12 q2m + 72 e - 12 w1 + 2
#   B23 = -45 p w1 + 15 p s1 - 3 q2 + 18 p
# Verified against the raw integrals of g_ab; then regular 2m-gon sums vs D17 targets; diamond m=2.
import sympy as sp

al, be, u, v, L, h, s0, s, c, m = sp.symbols('alpha beta u v L h s0 s c m', real=True)
w1 = al*u + be*v; e = al*be*u*v; q2m = al**2*u**2 + be**2*v**2
q2 = al**2 + be**2; p = al*be; s1 = al + be; r = h**2/L**2
A5  = 20*w1**3 - 20*e*w1 - 32*q2m - 48*e + 15*w1 - 2
B5  = -15*q2*w1 + 15*p*w1 + 5*p*s1 + 8*q2 - 12*p
A23 = -60*e*w1 + 12*q2m + 72*e - 12*w1 + 2
B23 = -45*p*w1 + 15*p*s1 - 3*q2 + 18*p
Z5law  = 360*al*be*h**2/L*(A5 + r*B5)
Z23law = 120*al*be*h**2/L*(A23 + r*B23)

# 1) match raw integrals from step_04
txt = open('steps/step_04_out.txt').read()
loc = {'alpha': al, 'beta': be, 'L': L, 'h': h, 's0': s0}
Z5raw  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals=loc)
Z23raw = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals=loc)
d5  = sp.expand(Z5raw  - Z5law.subs([(v, 1-u), (u, -s0/L)]))
d23 = sp.expand(Z23raw - Z23law.subs([(v, 1-u), (u, -s0/L)]))
print("law==raw  Z5:", sp.simplify(d5) == 0, " Z23:", sp.simplify(d23) == 0)

# 2) regular 2m-gon: alpha=beta=2 s^2, u=v=1/2, h_e = s, L = 2 s^2/c, h=pi/(2m), s=sin h, c=cos h
sub = [(al, 2*s**2), (be, 2*s**2), (u, sp.Rational(1,2)), (v, sp.Rational(1,2)),
       (h, s), (L, 2*s**2/c)]
def redc(expr):
    expr = sp.expand(expr)
    for pw in range(14, 1, -1):
        expr = sp.expand(expr.subs(c**pw, c**(pw-2)*(1-s**2)))
    return expr
Z5e  = redc(Z5law.subs(sub));  Z23e = redc(Z23law.subs(sub))
tgt5e  = 1440*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)   # = C_m^(5)/(2m)
tgt23e = -480*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)    # = C_m^(2,3)/(2m)
print("regular 2m-gon per-edge Z5 == target:",  redc(Z5e - tgt5e) == 0)
print("regular 2m-gon per-edge Z23 == target:", redc(Z23e - tgt23e) == 0)
# explicit m = 2..8: numeric check of the summed laws against D17 closed forms
import mpmath as mp
mp.mp.dps = 50
for mm in range(2, 9):
    hv = mp.pi/(2*mm); sv, cv = mp.sin(hv), mp.cos(hv)
    C5  = 2880*mm*sv**4*cv*(65*sv**6 - 63*sv**4 + 17*sv**2 - 1)
    C23 = -960*mm*cv*sv**4*(45*sv**6 - 27*sv**4 + 6*sv**2 - 1)
    # evaluate the edge laws at the regular data
    z5v  = sp.N(Z5law.subs([(al, 2*sv**2), (be, 2*sv**2), (u, sp.Rational(1,2)), (v, sp.Rational(1,2)),
                            (h, sv), (L, 2*sv**2/cv)]), 50)
    z23v = sp.N(Z23law.subs([(al, 2*sv**2), (be, 2*sv**2), (u, sp.Rational(1,2)), (v, sp.Rational(1,2)),
                             (h, sv), (L, 2*sv**2/cv)]), 50)
    print(f"m={mm}: |2m Z5e - C5| = {mp.nstr(abs(2*mm*mp.mpf(str(z5v)) - C5), 3)},"
          f" |2m Z23e - C23| = {mp.nstr(abs(2*mm*mp.mpf(str(z23v)) - C23), 3)}")
Cm5  = sp.expand(2*m*tgt5e);  Cm23 = sp.expand(2*m*tgt23e)
print("Sum_e Z5 = C_m^(5)  =", Cm5)
print("Sum_e Z23 = C_m^(2,3) =", Cm23)

# 3) diamond m=2: h=pi/4, s=c=sqrt2/2
sq = sp.sqrt(2)/2
val5  = sp.simplify(Cm5.subs([(s, sq), (c, sq), (m, 2)]))
val23 = sp.simplify(Cm23.subs([(s, sq), (c, sq), (m, 2)]))
print("diamond m=2: C_2^(5) =", val5, " (target -90 sqrt2:", sp.simplify(val5 + 90*sp.sqrt(2)) == 0, ")")
print("diamond m=2: C_2^(2,3) =", val23, " (target -210 sqrt2:", sp.simplify(val23 + 210*sp.sqrt(2)) == 0, ")")

with open('steps/step_05_out.txt', 'w') as f:
    f.write("Z_e^(5)  = 360 al be h^2/L [ A5 + (h^2/L^2) B5 ]\n")
    f.write("Z_e^(2,3)= 120 al be h^2/L [ A23 + (h^2/L^2) B23 ]\n")
    f.write(f"A5  = {A5}\nB5  = {B5}\nA23 = {A23}\nB23 = {B23}\n")
    f.write(f"law==raw Z5 {sp.simplify(d5)==0} Z23 {sp.simplify(d23)==0}\n")
    f.write(f"regular per-edge Z5==tgt {sp.expand(Z5e-tgt5e)==0}, Z23==tgt {sp.expand(Z23e-tgt23e)==0}\n")
    f.write(f"diamond C_2^(5)={val5}, C_2^(2,3)={val23}\n")
print("saved steps/step_05_out.txt")
