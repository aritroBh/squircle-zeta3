import sympy as sp
t, al, be, L, h, s0, u, v, r = sp.symbols('t alpha beta L h s0 u v r', real=True)
txt = open('steps/step_04_out.txt').read()
loc = {'alpha':al,'beta':be,'L':L,'h':h,'s0':s0,'t':t}
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals=loc)
Z23 = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals=loc)
swap = [(al, be), (be, al), (u, 1-u)]

def splitAB(Z):
    Ze  = sp.expand(Z.subs(s0, -u*L))
    Zn  = sp.expand(Ze/(al*be*h**2/L))
    A   = sp.expand(Zn.subs(h, 0))
    B   = sp.expand((Zn - A)*L**2/h**2)
    return A, B
A5, B5 = splitAB(Z5); A23, B23 = splitAB(Z23)
for nm, A, B in [("Z5", A5, B5), ("Z23", A23, B23)]:
    sA = sp.expand(A.subs(swap, simultaneous=True) - A)
    sB = sp.expand(B.subs(swap, simultaneous=True) - B)
    print(nm, "A symmetric:", sA==0, " B symmetric:", sB==0)

# Symmetric basis: w1 = al*u+be*v, w2 = al*u^2+be*v^2, w3 = al*u^3+be*v^3, e = al*be*u*v, s1=al+be, p=al*be, q2=al^2+be^2
w1 = al*u + be*v; w2 = al*u**2 + be*v**2; w3 = al*u**3 + be*v**3
e  = al*be*u*v
# note w3 = w1^3 - 3 e*w1  (since al^3 u^3+be^3 v^3 = (al u + be v)^3 - 3 al u be v (al u + be v))
# try to write A5/360 as combination: c1*w1^3 + c2*e*w1 + c3*w2 + c4*e + c5*w1 + c6
target = sp.expand(A5/360)
c1,c2,c3,c4,c5,c6 = sp.symbols('c1:7')
trial = c1*w1**3 + c2*e*w1 + c3*w2 + c4*e + c5*w1 + c6
sol = sp.solve(sp.Poly(sp.expand(trial - target).subs(v, 1-u), u, al, be).coeffs(), [c1,c2,c3,c4,c5,c6], dict=True)
print("A5/360 in {w1^3, e w1, w2, e, w1, 1}:", sol)
# B5: basis {al^3 u+be^3 v=w3b, al^2*be*v+al*be^2*u = e1b, al*be*(al+be)=p*s1, al^2+be^2=q2, al*be=p}
w3b = al**3*u + be**3*v; e1b = al**2*be*v + al*be**2*u
ps1 = al*be*(al+be); q2 = al**2+be**2; p = al*be
targetB = sp.expand(B5/360)
d1,d2,d3,d4,d5_,d6 = sp.symbols('d1:7')
trialB = d1*w3b + d2*e1b + d3*ps1 + d4*q2 + d5_*p + d6*(al*u**0)  # try w/o w1
solB = sp.solve(sp.Poly(sp.expand(trialB - targetB).subs(v, 1-u), u, al, be).coeffs(), [d1,d2,d3,d4,d5_,d6], dict=True)
print("B5/360 in {w3b, e1b, p s1, q2, p}:", solB)
