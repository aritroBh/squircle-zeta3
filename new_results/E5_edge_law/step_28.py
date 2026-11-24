import sympy as sp
al, be, u, v, L, h, s0 = sp.symbols('alpha beta u v L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
loc={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0}
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals=loc)
def splitAB(Z):
    Ze  = sp.expand(Z.subs(s0, -u*L))
    Zn  = sp.expand(Ze/(al*be*h**2/L))
    A   = sp.expand(Zn.subs(h, 0))
    B   = sp.expand((Zn - A)*L**2/h**2)
    return A, B
A5, B5 = splitAB(Z5)
print("A5 coeff al^3 u^3:", A5.coeff(al,3).coeff(u,3))
print("A5 coeff al^2 be u^3:", A5.coeff(al,2).coeff(be,1).coeff(u,3))
print("A5 coeff al u^3:", A5.coeff(al,1).coeff(be,0).coeff(u,3) if A5.coeff(al,1).coeff(u,3) else 0)
print(A5)
