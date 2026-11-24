import sympy as sp
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t})
P = {al: sp.Rational(3,7), be: sp.Rational(2,5), L: sp.Rational(4,3), h: sp.Rational(5,6), s0: sp.Rational(-2,9)}
lhs = Z5.subs(P)
Psw = {al: P[be], be: P[al], L: P[L], h: P[h], s0: -P[s0]-P[L]}
rhs = Z5.subs(Psw)
print("lhs - rhs =", sp.simplify(lhs - rhs))
D5 = sp.expand(Z5.subs({al: be, be: al, s0: -s0 - L}) - Z5)
print("D5(P) =", sp.simplify(D5.subs(P)))
# maybe D5 was built from a different Z5 string? check length
print(len(str(Z5)), len(str(D5)))
D5b = sp.expand(Z5.subs([(al, be), (be, al), (s0, -s0 - L)], simultaneous=True) - Z5)
print("D5b(P) =", sp.simplify(D5b.subs(P)), " D5-D5b:", sp.simplify(D5-D5b)==0)
