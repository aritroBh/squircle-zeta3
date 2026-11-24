
import pickle, os, json, sympy as sp, mpmath as mp
WORK = os.getcwd()
E3 = pickle.load(open(f"{WORK}/E3_fixed.pkl","rb"))          # FIXED pickle
Jtab = pickle.load(open(f"{WORK}/Jtab_t.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry.json"))
jloc_str = json.load(open(f"{WORK}/jloc.json"))
verif = {b: json.load(open(f"{WORK}/t3_{b}.json")) for b in ["b0_150","b1_110"]}

z2, z3, l2, l3 = sp.symbols('z2 z3 l2 l3')
Jloc = {tuple(map(int,k.split(","))): sp.sympify(v, locals={'z2':z2,'z3':z3}) for k,v in jloc_str.items()}

JL1,JL2,JL3,JAA,Jz2,Jx = sp.symbols('JL1 JL2 JL3 JAA Jz2 Jx')
ren = {sp.Symbol('L1'):JL1, sp.Symbol('L2'):JL2, sp.Symbol('L3'):JL3,
       sp.Symbol('AA'):JAA, sp.Symbol('z2'):Jz2, sp.Symbol('x'):Jx}
JtabR = {k: e.subs(ren, simultaneous=True) for k,e in Jtab.items()}

def logpq(rat):
    rat = sp.Rational(rat)
    a2 = sp.factorint(rat.p).get(2,0) - sp.factorint(rat.q).get(2,0)
    a3 = sp.factorint(rat.p).get(3,0) - sp.factorint(rat.q).get(3,0)
    assert sp.Rational(2)**a2 * sp.Rational(3)**a3 == rat, f"log of {rat} not in {{log2,log3}}"
    return a2*l2 + a3*l3

ATOMS = {
 '2':   {JL1:-l2, JL2:-z2/2, JL3:-sp.Rational(3,4)*z3, JAA:-sp.Rational(5,8)*z3, Jz2:z2},     # x=-1
 '1/2': {JL1: l2, JL2: z2/2 - l2**2/2, JL3: sp.Rational(7,8)*z3 - z2*l2/2 + l2**3/6,
         JAA: z3 - z2*l2/2, Jz2:z2},                                                          # x=1/2
}
symtab = {s.name: s for s in sp.symbols('A0 A1 B L2 V02 kap lA0 lA1 lB q0')}

def engine_arc(e):
    A0,A1,B = sp.Rational(e['A0']), sp.Rational(e['A1']), sp.Rational(e['B'])
    kap = sp.Rational(e['kap']); L2v = sp.Integer(e['L2']); q0 = sp.Integer(e['q0']); V02 = sp.Integer(e['V02'])
    w0,w1,wo = sp.Rational(e['w0']), sp.Rational(e['w1']), sp.Rational(e['wo'])
    sub = {symtab['A0']:A0, symtab['A1']:A1, symtab['B']:B, symtab['kap']:kap, symtab['L2']:L2v,
           symtab['q0']:q0, symtab['V02']:V02,
           symtab['lA0']:logpq(A0/w0), symtab['lA1']:logpq(A1/w1), symtab['lB']:logpq(B/wo)}  # UNWEIGHTED logs
    Np = sp.Integer(0)
    if kap != 0:
        for (n,a,b,c),co in E3['nl']:
            Np += co.subs(sub, simultaneous=True) * JtabR[(n,a,b,c)].subs({Jx:-kap})
        Np = sp.expand(Np)
    Lp = sp.Integer(0)
    for (n,a,b,c),co in E3['loc']:
        Lp += co.subs(sub, simultaneous=True) * Jloc[(n,a,b)]
    Lp = sp.expand(Lp)
    Earc_atoms = sp.expand(6*(Np + Lp))
    r = e['r']
    Earc_closed = sp.expand(Earc_atoms.subs(ATOMS[r], simultaneous=True)) if r in ATOMS else Earc_atoms
    assert sp.expand(Earc_closed - Earc_closed.subs({l3:0})) == 0, "log 3 did not cancel!"
    Earc_closed = Earc_closed.subs({l3:0})
    cz3 = sp.simplify(sp.expand(Earc_closed).coeff(z3))
    Ze = sp.sympify(e['Z'])
    Ne = sp.simplify(cz3 - Ze)
    return dict(Earc_closed=Earc_closed, cz3=cz3, Ze=Ze, Ne=Ne, L=sp.sqrt(L2v), kap=kap, r=r)

out = {}
mp.mp.dps = 120
for bkey, jname in [("0","b0_150"), ("1","b1_110")]:
    print("#"*92); print(f"# b = {bkey}   (engine: E3_fixed.pkl, unweighted logs)"); print("#"*92)
    arcs = verif[jname]['arcs']
    for e in geo[bkey]:
        res = engine_arc(e)
        print(f"--- b={bkey} edge {e['edge']}: r={res['r']}, kap={res['kap']}, atom x={-res['kap']} ---")
        print("E_arc closed form =", res['Earc_closed'])
        print("  = E/L basis: E = L*(", sp.expand(res['Earc_closed']/res['L']), ")")
        print("zeta(3) coeff cz3 =", res['cz3'])
        print("Z_e =", res['Ze'], "  =>  N_e = cz3 - Z_e =", res['Ne'])
        val = sp.N(res['Earc_closed'].subs({z3:sp.zeta(3), z2:sp.pi**2/6, l2:sp.log(2)}), 110)
        P3 = mp.mpf(arcs[e['edge']-1]['P3_arc'])
        diff = abs(mp.mpf(sp.N(val,100)) - P3)
        print("engine E_arc =", mp.nstr(mp.mpf(sp.N(val,70)), 70))
        print("verifier P3  =", mp.nstr(P3, 70))
        print("abs diff =", mp.nstr(diff, 3))
        rel = arcs[e['edge']-1]['relation']
        X = res['Earc_closed']/res['L']
        basis = [X, z3, 6*z2*l2, l2**3, 6*z2, l2**2, l2, sp.Integer(1)]
        exactres = sp.expand(sum(sp.Integer(c)*bb for c,bb in zip(rel,basis)))
        print("verifier relation", rel, "-> exact residual:", exactres)
        out[f"b{bkey}_e{e['edge']}"] = dict(Earc_closed=str(res['Earc_closed']), cz3=str(res['cz3']),
                                            Ze=str(res['Ze']), Ne=str(res['Ne']), r=str(res['r']),
                                            numeric=mp.nstr(mp.mpf(sp.N(val,80)),80), diff=mp.nstr(diff,3),
                                            rel_residual=str(exactres))
json.dump(out, open(f"{WORK}/engine_arcs_fixed.json","w"), indent=1)
print("saved engine_arcs_fixed.json")
