import sympy as sp, itertools
al, be, u, v, L, h, s0 = sp.symbols('alpha beta u v L h s0', real=True)
txt = open('steps/step_04_out.txt').read()
loc={'alpha':al,'beta':be,'L':L,'h':h,'s0':s0}
Z5  = sp.sympify(txt.split('Z5 raw = ')[1].split('\n')[0], locals=loc)
Z23 = sp.sympify(txt.split('Z23 raw = ')[1].split('\n')[0], locals=loc)
def splitAB(Z):
    Ze  = sp.expand(Z.subs(s0, -u*L)); Zn  = sp.expand(Ze/(al*be*h**2/L))
    return sp.expand(Zn.subs(h, 0)), sp.expand((Zn - sp.expand(Zn.subs(h,0)))*L**2/h**2)
A5, B5 = splitAB(Z5); A23, B23 = splitAB(Z23)
swap = [(al, be), (be, al), (u, v)]
# build library of symmetric polynomials: all symmetrised monomials al^a be^b u^c (v^d = (v)) with a+b<=3, c<=3
lib = {}
seen = set()
for a in range(4):
    for b in range(4-a):
        for c in range(4):
            for dd in range(4):
                if (a,b,c,dd) in seen: continue
                m = al**a*be**b*u**c*v**dd
                ms = m + m.subs(swap, simultaneous=True)
                key = sp.expand(ms.subs(v,1-u))
                if key == 0: continue
                lib[(a,b,c,dd)] = key
                seen.add((a,b,c,dd)); seen.add((b,a,dd,c))
names = list(lib); Bm = sp.Matrix([ [sp.Poly(lib[n], u, al, be).coeff_monomial(mn) for mn in
        sorted(set().union(*[set(sp.Poly(lib[k],u,al,be).monoms()) for k in names]))] for n in names])
monoms = sorted(set().union(*[set(sp.Poly(lib[k],u,al,be).monoms()) for k in names]))
Bm = sp.Matrix([[sp.Poly(lib[n], u, al, be).coeff_monomial(mn) for mn in monoms] for n in names])
def fitS(expr, scale):
    pe = sp.Poly(sp.expand(expr/scale).subs(v,1-u) if not sp.Poly(expr,u,al,be).terms() else expr/scale, u, al, be)
    pe = sp.Poly(sp.expand(sp.expand(expr).subs(v,1-u)/scale), u, al, be)
    vec = sp.Matrix([pe.coeff_monomial(mn) for mn in monoms])
    sol = sp.linsolve((Bm.T, vec))
    return dict(zip(names, list(list(sol)[0])))
for nm, A, B, sc in [("Z5", A5, B5, 360), ("Z23", A23, B23, 120)]:
    fA = fitS(A, sc); fB = fitS(B, sc)
    fA = {k:v_ for k,v_ in fA.items() if v_ != 0}; fB = {k:v_ for k,v_ in fB.items() if v_ != 0}
    print(nm, "A/%d ="%sc, fA)
    print(nm, "B/%d ="%sc, fB)
