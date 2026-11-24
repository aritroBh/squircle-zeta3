import sympy as sp
al, be, u, v, L, h = sp.symbols('alpha beta u v L h', real=True)
txt = open('steps/step_04_out.txt').read()
s0 = sp.symbols('s0', real=True)
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0})
Z23 = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0})
w1 = al*u + be*v; w2 = al*u**2 + be*v**2; e = al*be*u*v
r = h**2/L**2
# candidate closed forms
A5f = 20*w1**3 - 20*e*w1 - 32*w2 - 48*e + 15*w1 - 2
B5f = -15*(al**3*u + be**3*v) - 30*(al**2*be*v + al*be**2*u) + 20*al*be*(al+be) + 8*(al**2+be**2) - 12*al*be
Z5f = 360*al*be*h**2/L*(A5f + r*B5f)
A23f = sp.expand((sp.expand(Z23.subs(s0,-u*L))/(al*be*h**2/L)).subs(h,0))
B23f = sp.expand((sp.expand(Z23.subs(s0,-u*L))/(al*be*h**2/L) - A23f)*L**2/h**2)
print("Z5 candidate match:", sp.expand(Z5.subs(s0, -u*L) - Z5f.subs(v, 1-u)) == 0)
print("A23/120 =", sp.factor(A23f/120))
print("B23/120 =", sp.factor(B23f/120))
