import sympy as sp, pickle, random
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t})
Z23 = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t})
random.seed(1)
for _ in range(3):
    va, vb, vL, vh, vs = [sp.Rational(random.randint(1,9), random.randint(1,9)) for _ in range(5)]
    d5  = sp.simplify(Z5.subs({al:va,be:vb,L:vL,h:vh,s0:vs})  - Z5.subs({al:vb,be:va,L:vL,h:vh,s0:-vs-vL}))
    d23 = sp.simplify(Z23.subs({al:va,be:vb,L:vL,h:vh,s0:vs}) - Z23.subs({al:vb,be:va,L:vL,h:vh,s0:-vs-vL}))
    print("raw reversal-invariance:", d5 == 0, d23 == 0)
# symbolic reversal check on Z5
D5 = sp.expand(Z5.subs({al: be, be: al, s0: -s0 - L}) - Z5)
print("Z5 symbolic reversal defect:", sp.factor(D5))
D23 = sp.expand(Z23.subs({al: be, be: al, s0: -s0 - L}) - Z23)
print("Z23 symbolic reversal defect:", sp.factor(D23))
