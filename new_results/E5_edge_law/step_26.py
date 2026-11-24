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

# symmetric basis elements (u only; v -> 1-u inside)
def sym_basis():
    m1 = al*u + be*v; m2 = al*u**2 + be*v**2; m3 = al*u**3 + be*v**3
    mm = al*be*u*v; em1 = al*be*u*v*(al*u + be*v)
    n31 = al**3*u + be**3*v; n22 = al**2*be*v + al*be**2*u
    ps1 = al*be*(al+be); q2 = al**2+be**2; p = al*be
    n21 = al**2*u + be**2*v      # deg (2,1)
    m2m = al*be*(u+v)            # = al*be on locus
    return dict(m1=m1, m2=m2, m3=m3, mm=mm, em1=em1, n31=n31, n22=n22, ps1=ps1, q2=q2, p=p, n21=n21)
SB = sym_basis()
def fit(expr, names):
    cs = sp.symbols('c:'+str(len(names)))
    trial = sum(c*SB[n] for c, n in zip(cs, names))
    diff = sp.Poly(sp.expand((trial - expr).subs(v, 1-u)), u, al, be)
    sol = sp.solve(diff.coeffs(), cs, dict=True)
    return sol
print("A5/360 fit {m3,em1,m2,mm,m1,1}:", fit(A5/360, ['m3','em1','m2','mm','m1']))
# include constant: '1'
names5 = ['m3','em1','m2','mm','m1']
cs = sp.symbols('c:6'); SBv = [SB[n] for n in names5]+[sp.Integer(1)]
diff = sp.Poly(sp.expand((sum(c*bb for c,bb in zip(cs,SBv)) - A5/360).subs(v,1-u)), u, al, be)
sol = sp.solve(diff.coeffs(), cs, dict=True); print("A5/360 with const:", sol)
cs = sp.symbols('c:6'); SBv = [SB[n] for n in ['n31','n22','ps1','q2','p']]+[sp.Integer(1)]
diff = sp.Poly(sp.expand((sum(c*bb for c,bb in zip(cs,SBv)) - B5/360).subs(v,1-u)), u, al, be)
print("B5/360:", sp.solve(diff.coeffs(), cs, dict=True))
cs = sp.symbols('c:7'); SBv = [SB[n] for n in ['em1','m2','mm','m1','n21','p']]+[sp.Integer(1)]
diff = sp.Poly(sp.expand((sum(c*bb for c,bb in zip(cs,SBv)) - A23/120).subs(v,1-u)), u, al, be)
print("A23/120:", sp.solve(diff.coeffs(), cs, dict=True))
cs = sp.symbols('c:5'); SBv = [SB[n] for n in ['n21','p','q2']]+[sp.Integer(1)]
diff = sp.Poly(sp.expand((sum(c*bb for c,bb in zip(cs,SBv)) - B23/120).subs(v,1-u)), u, al, be)
print("B23/120:", sp.solve(diff.coeffs(), cs, dict=True))
