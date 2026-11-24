
import pickle, json, sympy as sp, time, os
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
E3 = pickle.load(open(f"{WORK}/E3_fixed.pkl","rb"))
Jtab = pickle.load(open(f"{WORK}/Jtab_t.pkl","rb"))
jloc = json.load(open(f"{WORK}/jloc.json"))

z2,z3 = sp.symbols('z2 z3')
Jloc = {tuple(map(int,k.split(","))): sp.sympify(v, locals={'z2':z2,'z3':z3}) for k,v in jloc.items()}

A0,A1,B,kap,L2s,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
L1,L2,L3,L3L,AA,x = sp.symbols('L1 L2 L3 L3L AA x')
ren = {sp.Symbol('L1'):L1, sp.Symbol('L2'):L2, sp.Symbol('L3'):L3, sp.Symbol('AA'):AA, sp.Symbol('z2'):z2, sp.Symbol('x'):x}
JtabR = {k: sp.expand(e.subs(ren, simultaneous=True)) for k,e in Jtab.items()}

pars = {'A0':A0,'A1':A1,'B':B,'kap':kap,'L2':L2s,'q0':q0,'V02':V02,'lA0':lA0,'lA1':lA1,'lB':lB}
symtab = {s: sp.Symbol(s) for s in pars}

# non-local part: substitute x -> -kap in J-table
Np = sp.Integer(0)
for (n,a,b,c),co in E3['nl']:
    co2 = co.subs({symtab[k]:v for k,v in pars.items()}, simultaneous=True)
    Np += co2 * JtabR[(n,a,b,c)].subs({x:-kap}, simultaneous=True)
Np = sp.expand(Np)
print("Np terms:", len(sp.Add.make_args(Np)), f"({time.time()-t0:.0f}s)")
Lp = sp.Integer(0)
for (n,a,b,c),co in E3['loc']:
    co2 = co.subs({symtab[k]:v for k,v in pars.items()}, simultaneous=True)
    Lp += co2 * Jloc[(n,a,b)]
Lp = sp.expand(Lp)
print("Lp terms:", len(sp.Add.make_args(Lp)), f"({time.time()-t0:.0f}s)")
Earc = sp.expand(6*(Np+Lp))
# eliminate AA:  A(x) = 2 Li3(x) + Li3(x/(x-1)) + Li2(x)*L1 + L1^3/6   (ln(1-x) = -L1)
Earc = sp.expand(Earc.subs({AA: 2*L3 + L3L + L2*L1 + L1**3/sp.Integer(6)}))
print("Earc terms after AA elimination:", len(sp.Add.make_args(Earc)), f"({time.time()-t0:.0f}s)")
poly = sp.Poly(Earc, L1,L2,L3,L3L,z2,z3)
mons = {}
for exps, coef in poly.terms():
    mons[exps] = sp.expand(coef)
print("distinct monomials:", len(mons))
namemap = {(0,0,1,0,0,0):'Li3(x)', (0,0,0,1,0,0):'Li3(x/(x-1))', (1,1,0,0,0,0):'Li2(x)*ln(1-x)',
           (3,0,0,0,0,0):'ln^3(1-x)', (0,1,0,0,0,0):'Li2(x)', (2,0,0,0,0,0):'ln^2(1-x)',
           (1,0,0,0,0,0):'ln(1-x)', (0,0,0,0,0,1):'zeta3', (1,0,0,0,1,0):'zeta2*ln(1-x)',
           (0,0,0,0,1,0):'zeta2', (0,0,0,0,0,0):'1'}
extra = [e for e in mons if e not in namemap]
print("monomials outside the task basis:", extra)
out = {namemap[e]: str(c) for e,c in mons.items() if e in namemap}
pickle.dump({namemap.get(e, str(e)): c for e,c in mons.items()}, open(f"{WORK}/Earc_basis_raw.pkl","wb"))
json.dump(out, open(f"{WORK}/Earc_basis_raw.json","w"), indent=1)
for e,c in sorted(mons.items(), key=lambda kv: str(kv[0])):
    print(namemap.get(e,str(e)), "coef #terms:", len(sp.Add.make_args(c)))
print("saved Earc_basis_raw.pkl", f"({time.time()-t0:.0f}s)")
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_07: built E_arc = 6(Np+Lp) symbolically, eliminated AA via 6A=12Li3+6Li3L-6Li2 ln-ln^3; collected on basis; {len(mons)} monomials, extras={extra}. Saved Earc_basis_raw.pkl.\n")
