
import pickle, json, sympy as sp, time, os
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
E3 = pickle.load(open(f"{WORK}/E3_fixed.pkl","rb"))
Jtab = pickle.load(open(f"{WORK}/Jtab_t.pkl","rb"))
jloc = json.load(open(f"{WORK}/jloc.json"))

z2,z3 = sp.symbols('z2 z3')
Jloc = {tuple(map(int,k.split(","))): sp.sympify(v, locals={'z2':z2,'z3':z3}) for k,v in jloc.items()}

# parameters
A0,A1,B,kap,Lp2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# atoms: l1=-ln(1-x), l2a=Li2(x), l3a=Li3(x), l3l=Li3(x/(x-1)), AA=A(x)
l1,l2a,l3a,l3l,AA,x = sp.symbols('l1 l2a l3a l3l AA x')
ren = {sp.Symbol('L1'):l1, sp.Symbol('L2'):l2a, sp.Symbol('L3'):l3a, sp.Symbol('AA'):AA, sp.Symbol('z2'):z2, sp.Symbol('x'):x}
JtabR = {k: sp.expand(e.subs(ren, simultaneous=True)) for k,e in Jtab.items()}

symtab = {s: sp.Symbol(s) for s in ['A0','A1','B','kap','L2','q0','V02','lA0','lA1','lB']}
par_sub = {symtab['A0']:A0, symtab['A1']:A1, symtab['B']:B, symtab['kap']:kap, symtab['L2']:Lp2,
           symtab['q0']:q0, symtab['V02']:V02, symtab['lA0']:lA0, symtab['lA1']:lA1, symtab['lB']:lB}

Np = sp.Integer(0)
for (n,a,b,c),co in E3['nl']:
    Np += co.subs(par_sub, simultaneous=True) * JtabR[(n,a,b,c)].subs({x:-kap}, simultaneous=True)
Np = sp.expand(Np)
Lp = sp.Integer(0)
for (n,a,b,c),co in E3['loc']:
    Lp += co.subs(par_sub, simultaneous=True) * Jloc[(n,a,b)]
Lp = sp.expand(Lp)
Earc = sp.expand(6*(Np+Lp))
Earc = sp.expand(Earc.subs({AA: 2*l3a + l3l + l2a*l1 + l1**3/sp.Integer(6)}))
print("Earc terms:", len(sp.Add.make_args(Earc)), f"({time.time()-t0:.0f}s)")
poly = sp.Poly(Earc, l1,l2a,l3a,l3l,z2,z3)
mons = {exps: sp.expand(coef) for exps, coef in poly.terms()}
namemap = {(0,0,1,0,0,0):'Li3(x)', (0,0,0,1,0,0):'Li3(x/(x-1))', (1,1,0,0,0,0):'Li2(x)*ln(1-x)',
           (3,0,0,0,0,0):'ln^3(1-x)', (0,1,0,0,0,0):'Li2(x)', (2,0,0,0,0,0):'ln^2(1-x)',
           (1,0,0,0,0,0):'ln(1-x)', (0,0,0,0,0,1):'zeta3', (1,0,0,0,1,0):'zeta2*ln(1-x)',
           (0,0,0,0,1,0):'zeta2', (0,0,0,0,0,0):'1'}
extra = [e for e in mons if e not in namemap]
print("monomials outside the task basis:", extra)
pickle.dump({namemap.get(e, str(e)): c for e,c in mons.items()}, open(f"{WORK}/Earc_basis_raw.pkl","wb"))
for e,c in sorted(mons.items(), key=lambda kv: str(kv[0])):
    print(f"{namemap.get(e,str(e)):>16}: coef #terms {len(sp.Add.make_args(c))}")
print(f"saved ({time.time()-t0:.0f}s)")
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_08: built E_arc=6(Np+Lp) symbolically (1613 terms), eliminated AA via 6A=12Li3+6Li3L-6Li2*ln-ln^3; collected on 11-element basis; extras outside basis: {extra}. Saved Earc_basis_raw.pkl.\n")
