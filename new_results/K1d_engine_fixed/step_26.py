import pickle, json, sympy as sp
E3 = pickle.load(open('E3_fixed.pkl','rb'))
Jtab = pickle.load(open('Jtab_t.pkl','rb'))
geo = json.load(open('geometry.json'))
jloc_str = json.load(open('jloc.json'))
z2,z3,l2,l3 = sp.symbols('z2 z3 l2 l3')
Jloc = {tuple(map(int,k.split(","))): sp.sympify(v, locals={'z2':z2,'z3':z3}) for k,v in jloc_str.items()}
JL1,JL2,JL3,JAA,Jz2,Jx = sp.symbols('JL1 JL2 JL3 JAA Jz2 Jx')
ren = {sp.Symbol('L1'):JL1, sp.Symbol('L2'):JL2, sp.Symbol('L3'):JL3, sp.Symbol('AA'):JAA, sp.Symbol('z2'):Jz2, sp.Symbol('x'):Jx}
JtabR = {k: e.subs(ren, simultaneous=True) for k,e in Jtab.items()}
def logpq(rat):
    rat = sp.Rational(rat)
    a2 = sp.factorint(rat.p).get(2,0) - sp.factorint(rat.q).get(2,0)
    a3 = sp.factorint(rat.p).get(3,0) - sp.factorint(rat.q).get(3,0)
    return a2*l2 + a3*l3
symtab = {s.name: s for s in sp.symbols('A0 A1 B L2 V02 kap lA0 lA1 lB q0')}

for bkey in ['0','1']:
    for e in geo[bkey]:
        A0,A1,B = sp.Rational(e['A0']), sp.Rational(e['A1']), sp.Rational(e['B'])
        kap = sp.Rational(e['kap']); L2v = sp.Integer(e['L2']); q0 = sp.Integer(e['q0']); V02 = sp.Integer(e['V02'])
        w0,w1,wo = sp.Rational(e['w0']), sp.Rational(e['w1']), sp.Rational(e['wo'])
        sub = {symtab['A0']:A0, symtab['A1']:A1, symtab['B']:B, symtab['kap']:kap, symtab['L2']:L2v,
               symtab['q0']:q0, symtab['V02']:V02,
               symtab['lA0']:logpq(A0/w0), symtab['lA1']:logpq(A1/w1), symtab['lB']:logpq(B/wo)}
        # atom coefficients of 6*Np (x kept symbolic)
        PhiA = sp.Integer(0); PhiL = sp.Integer(0)
        for (n,a,b,c),co in E3['nl']:
            Je = JtabR[(n,a,b,c)]
            cA = sp.expand(Je).coeff(JAA); cL = sp.expand(Je).coeff(JL3)
            cs = co.subs(sub, simultaneous=True)
            PhiA += cs*cA; PhiL += cs*cL
        PhiA = sp.simplify(sp.expand(PhiA)*6); PhiL = sp.simplify(sp.expand(PhiL)*6)
        # local z3
        Lp = sp.Integer(0)
        for (n,a,b,c),co in E3['loc']:
            Lp += co.subs(sub, simultaneous=True) * Jloc[(n,a,b)]
        z3loc = sp.simplify(sp.expand(6*Lp).coeff(z3))
        Ze = sp.sympify(e['Z'])
        print(f"b={bkey} edge {e['edge']} r={e['r']} kap={e['kap']}:")
        print("   Phi_A(x) =", PhiA)
        print("   Phi_L(x) =", PhiL)
        print("   at x=-kap:", sp.simplify(PhiA.subs(Jx,-kap)), " | ", sp.simplify(PhiL.subs(Jx,-kap)))
        print("   z3-coeff 6*Lp =", z3loc, "  Z_e =", Ze, "  equal:", sp.simplify(z3loc-Ze)==0)
