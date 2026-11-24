import sympy as sp, random
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t})
random.seed(1)
for _ in range(3):
    va, vb, vL, vh, vs = [sp.Rational(random.randint(1,9), random.randint(1,9)) for _ in range(5)]
    lhs = Z5.subs({al:va,be:vb,L:vL,h:vh,s0:vs})
    rhs = Z5.subs({al:vb,be:va,L:vL,h:vh,s0:-vs-vL})
    print("vals:", va, vb, vL, vh, vs, " lhs:", lhs, " rhs:", rhs, " equal:", sp.simplify(lhs-rhs)==0)
