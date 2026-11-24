import sympy as sp
al, be, u, v, L, h, s0 = sp.symbols('alpha beta u v L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
z23str = txt.split('Z23 raw = ')[1].split('\n')[0]
print("Z23 raw head:", z23str[:300])
Z23 = sp.sympify(z23str, locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0})
Zn = sp.together(sp.expand(Z23.subs(s0, -u*L))/(al*be*h**2/L))
print("Zn =", sp.factor(Zn))
