import sympy as sp
al, be, u, v, L, h, s0 = sp.symbols('alpha beta u v L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
loc={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0}
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals=loc)
Z23 = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals=loc)
def splitAB(Z):
    Ze  = sp.expand(Z.subs(s0, -u*L))
    Zn  = sp.expand(Ze/(al*be*h**2/L))
    A   = sp.expand(Zn.subs(h, 0))
    B   = sp.expand((Zn - A)*L**2/h**2)
    return A, B
A5, B5 = splitAB(Z5); A23, B23 = splitAB(Z23)
# my hand guess
m1 = al*u + be*v; m2 = al*u**2 + be*v**2; m3 = al*u**3 + be*v**3
mm = al*be*u*v; em1 = mm*m1
guess = 360*(20*m3 + 40*em1 - 32*m2 - 48*mm + 15*m1 - 2)
print("A5 - guess:", sp.expand((A5 - guess).subs(v, 1-u)))
print("A23/120 =", sp.collect(A23/120, u))
print("B23/120 =", sp.factor(B23/120))
