import sympy as sp, pickle
t, al, be, L, h, s0, u, v, r = sp.symbols('t alpha beta L h s0 u v r', real=True)
Z5u  = sp.sympify(open('steps/step_04_out.txt').read().split('Z5 raw = ')[1].split('\n')[0],
                  locals={'t':t,'alpha':al,'beta':be,'L':L,'h':h,'s0':s0})
Z23u = sp.sympify(open('steps/step_04_out.txt').read().split('Z23 raw = ')[1].split('\n')[0],
                  locals={'t':t,'alpha':al,'beta':be,'L':L,'h':h,'s0':s0})
# substitute s0 = -uL; then express in u only, separating h^2/L^2 = r parts
Z5e  = sp.expand(Z5u.subs(s0, -u*L))
Z23e = sp.expand(Z23u.subs(s0, -u*L))
A5  = sp.expand(Z5e.coeff(h, 2).subs(L, 1)/360*L**0)  # placeholder
# Extract: Z = al*be*h^2/L * [polyA + r*polyB], r = h^2/L^2
def splitAB(Z):
    Zn = sp.expand(Z/(al*be*h**2/L))
    P  = sp.Poly(sp.expand(Zn*L**2), L)          # degree in L
    A  = sp.expand(Zn.subs(h, 0))                # part with no extra h
    B  = sp.expand((Zn - A)*L**2/h**2)
    return sp.expand(A), sp.expand(B)
A5, B5 = splitAB(Z5e); A23, B23 = splitAB(Z23e)
print("A5/360 =", sp.factor(A5/360))
print("B5/360 =", sp.factor(B5/360))
print("A23/120 =", sp.factor(A23/120))
print("B23/120 =", sp.factor(B23/120))
# symmetry check under (al,u) <-> (be,1-u)
for nm, A, B in [("Z5", A5, B5), ("Z23", A23, B23)]:
    symA = sp.expand(A.subs({al: be, be: al, u: 1-u}) - A)
    symB = sp.expand(B.subs({al: be, be: al, u: 1-u}) - B)
    print(nm, "symmetric A:", symA == 0, " symmetric B:", symB == 0)
