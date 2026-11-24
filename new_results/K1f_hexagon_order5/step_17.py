
import sympy as sp, pickle
WD = '<project>/research/08_theory/K1f_hexagon_order5'
K1D = '<project>/research/08_theory/K1d_engine_fixed'
A0,A1,B,kap = sp.symbols('A0 A1 B kap')
L2,q0,V02 = sp.symbols('L2 q0 V02')
lA0,lA1,lB = sp.symbols('lA0 lA1 lB')
t,X,Y,W = sp.symbols('t X Y W')
local = {'A0':A0,'A1':A1,'B':B,'kap':kap,'L2':L2,'q0':q0,'V02':V02,'lA0':lA0,'lA1':lA1,'lB':lB,'t':t,'X':X,'Y':Y,'W':W}
Ea = pickle.load(open(f'{WD}/E5_general_approachA.pkl','rb'))['E']
E3 = sp.sympify(Ea[3], locals=local)
P = sp.Poly(sp.expand(E3), t, X, Y, W)
my = {tuple(m): sp.expand(c) for m,c in P.terms()}
e3 = pickle.load(open(f'{K1D}/E3_fixed.pkl','rb'))
k1d = {}
for mon,co in e3['nl']+e3['loc']:
    k1d[tuple(mon)] = sp.expand(co)
# compare
allk = set(my)|set(k1d)
bad=0
for k in sorted(allk):
    diff = sp.simplify(my.get(k,sp.Integer(0)) - k1d.get(k,sp.Integer(0)))
    if diff!=0:
        bad+=1
        if bad<=10: print("MISMATCH", k, "mine=",my.get(k), "k1d=",k1d.get(k))
print("total monomials mine:", len(my), "k1d:", len(k1d), "mismatches:", bad)
