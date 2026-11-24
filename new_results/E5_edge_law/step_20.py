import sympy as sp
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t})
D5 = sp.expand(Z5.subs({al: be, be: al, s0: -s0 - L}) - Z5)
val = sp.simplify(D5.subs({al: sp.Rational(3,7), be: sp.Rational(2,5), L: sp.Rational(4,3), h: sp.Rational(5,6), s0: sp.Rational(-2,9)}))
print("defect at random point:", val)
print("defect == 0 after simplify:", sp.simplify(D5) == 0)
