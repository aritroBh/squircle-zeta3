code = r'''# step_05.py: THE ZETA(5) EDGE LAW -- closed forms and symbolic verification.
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
Z5e  = sp.expand(Z5law.subs(sub));  Z23e = sp.expand(Z23law.subs(sub))
Z5e  = sp.expand(Z5e.subs(c**2, 1 - s**2));  Z23e = sp.expand(Z23e.subs(c**2, 1 - s**2))
tgt5e  = 1440*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)   # = C_m^(5)/(2m)
tgt23e = -480*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)    # = C_m^(2,3)/(2m)
print("regular 2m-gon per-edge Z5 == target:",  sp.expand(Z5e - tgt5e) == 0)
print("regular 2m-gon per-edge Z23 == target:", sp.expand(Z23e - tgt23e) == 0)
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
'''
open('steps/step_05.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_05.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
