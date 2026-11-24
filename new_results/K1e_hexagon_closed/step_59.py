
import pickle, sympy as sp, mpmath as mp, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps=80
t0=time.time()
E3 = pickle.load(open(f"{WORK}/E3_fixed.pkl","rb"))
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
L2s=sp.Symbol('L2')
pars_syms = sp.symbols('A1 B kap q0 V02 lA0 lA1 lB')
atom_syms = sp.symbols('l3a l3l l2l1 l1_3 l2a l1_2 l1 z3 z2l1 z2 one')
basis_order = ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','ln^3(1-x)','Li2(x)','ln^2(1-x)','ln(1-x)','zeta3','zeta2*ln(1-x)','zeta2','1']
Eexpr = sum(num/L2s**sp.Rational(5,2)*atom for (num,den),atom in zip([clean[n] for n in basis_order], atom_syms))
f = sp.lambdify(pars_syms + (L2s,) + atom_syms, Eexpr, 'mpmath')
# pickle integrand as a python function of t given numeric params
symtab = {s: sp.Symbol(s) for s in ['A0','A1','B','kap','L2','q0','V02','lA0','lA1','lB']}
z2v=mp.pi**2/6; z3v=mp.zeta(3)
tests = [
    dict(A1=sp.Rational(2,3), B=sp.Rational(1,3), kap=sp.Rational(4,3), L2=sp.Integer(29), q0=sp.Integer(-7), V02=sp.Integer(10), la=(2,3,5)),   # r=7/3, x=-4/3
    dict(A1=sp.Rational(3,5), B=sp.Rational(2,5), kap=sp.Rational(-3,4), L2=sp.Integer(17), q0=sp.Integer(3), V02=sp.Integer(5), la=(5,2,3)),   # r=1/4, x=3/4
    dict(A1=sp.Rational(1,2), B=sp.Rational(1,2), kap=sp.Rational(3), L2=sp.Integer(41), q0=sp.Integer(-11), V02=sp.Integer(13), la=(3,5,2)),   # r=4, x=-3
]
for i,tc in enumerate(tests):
    A1,B,kap,L2v,q0v,V02v = tc['A1'],tc['B'],tc['kap'],tc['L2'],tc['q0'],tc['V02']
    A0 = A1 - B*kap
    assert A1+B==1 and L2v*V02v-q0v**2>0
    la0,la1,lb = [mp.log(p) for p in tc['la']]
    subs = {symtab['A0']:A0, symtab['A1']:A1, symtab['B']:B, symtab['kap']:kap, symtab['L2']:L2v,
            symtab['q0']:q0v, symtab['V02']:V02v, symtab['lA0']:0, symtab['lA1']:0, symtab['lB']:0}
    # build numeric integrand (logs handled separately)
    pieces=[]
    for mono,co in E3['nl']+E3['loc']:
        n,a,b,c = mono
        co_num = co.subs({k:v for k,v in subs.items()})
        # substitute logs numerically via N with lA symbols
        co_n = sp.N(co_num.subs({symtab['lA0']:la0, symtab['lA1']:la1, symtab['lB']:lb}), 60)
        pieces.append((n,a,b,c, mp.mpf(str(co_n))))
    kapn=mp.mpf(kap)
    def integrand(t):
        if t==0: t=mp.mpf('1e-300')
        if t==1: t=1-mp.mpf('1e-300')
        lt=mp.log(t); l1t=mp.log(1-t); lkt=mp.log(1+kapn*t)
        s=mp.mpf(0)
        for n,a,b,c,cn in pieces:
            s += cn * t**n * lt**a * l1t**b * lkt**c
        return s
    # quadrature with endpoint splitting
    I = mp.quad(integrand, [0, mp.mpf('0.5'), 1], method='tanh-sinh')
    Eq = 6*I
    # closed form
    x = -mp.mpf(kap)
    l1=-mp.log(1-x); l2a=mp.polylog(2,x); l3a=mp.polylog(3,x); l3l=mp.polylog(3,x/(x-1))
    val = f(mp.mpf(A1),mp.mpf(B),mp.mpf(kap),mp.mpf(q0v),mp.mpf(V02v),la0,la1,lb,mp.mpf(L2v),
            l3a,l3l,l2a*l1,l1**3,l2a,l1**2,l1,z3v,z2v*l1,z2v,mp.mpf(1))
    d=abs(val-Eq)
    print(f"test {i+1} (r={kap+1}, x={mp.nstr(x,6)}): quad E_arc={mp.nstr(Eq,25)}")
    print(f"          closed form   ={mp.nstr(val,25)}   |diff|={mp.nstr(d,3)}  ({time.time()-t0:.0f}s)")
