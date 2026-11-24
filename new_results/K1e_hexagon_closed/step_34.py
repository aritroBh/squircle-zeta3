
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
# (a) every non-local J entry vanishes at x=0 (series check)
Jtab = pickle.load(open(f"{WORK}/Jtab_t.pkl","rb"))
l1,l2a,l3a,AA,z2,x = sp.symbols('l1 l2a l3a AA z2 x')
ren = {sp.Symbol('L1'):l1, sp.Symbol('L2'):l2a, sp.Symbol('L3'):l3a, sp.Symbol('AA'):AA, sp.Symbol('z2'):z2, sp.Symbol('x'):x}
# series of atoms at x=0 to O(x^4)
aser = {l1: x+x**2/2+x**3/3+x**4/4, l2a: x+x**2/4+x**3/9+x**4/16, l3a: x+x**2/8+x**3/27+x**4/64,
        AA: x+3*x**2/8+11*x**3/54+25*x**4/192}
bad=[]
for k,e in Jtab.items():
    se = sp.expand(e.subs(ren, simultaneous=True))
    se = sp.expand(se.subs(aser))
    se = sp.series(se, x, 0, 1).removeO().expand()
    const = sp.expand(se).subs({x:0})
    if sp.simplify(const)!=0: bad.append((k,const))
print("(a) J entries with nonzero constant term at x=0:", bad)
# (b) E_clean at kap=0 with atoms->0 equals 6*Lp at kap=0
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
E3 = pickle.load(open(f"{WORK}/E3_fixed.pkl","rb"))
import json
jloc = json.load(open(f"{WORK}/jloc.json"))
Jloc = {tuple(map(int,k.split(","))): sp.sympify(v, locals={'z2':sp.Symbol('z2'),'z3':sp.Symbol('z3')}) for k,v in jloc.items()}
A1,B,kap,q0,V02,L2,lA0,lA1,lB,z3 = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB z3')
symtab = {s: sp.Symbol(s) for s in ['A0','A1','B','kap','L2','q0','V02','lA0','lA1','lB']}
Lp6 = sp.Integer(0)
for (n,a,b,c),co in E3['loc']:
    Lp6 += co*Jloc[(n,a,b)]
Lp6 = sp.expand(6*Lp6)
Lp6 = sp.expand(Lp6.subs({symtab['A0']:A1, symtab['kap']:0}))  # kap=0 => A0=A1
E0 = sp.Integer(0)
for name in ['1','zeta2','zeta3']:
    num,den = clean[name]
    E0 += sp.expand(num.subs({kap:0, sp.Symbol('A0'):A1}))/L2**sp.Rational(5,2)*{'1':1,'zeta2':sp.Symbol('z2'),'zeta3':z3}[name]
E0 = sp.expand(E0)
Lp6n = sp.expand(Lp6.subs({symtab['L2']:L2, symtab['q0']:q0, symtab['V02']:V02, symtab['lA0']:lA0, symtab['lA1']:lA1, symtab['lB']:lB, symtab['A1']:A1, symtab['B']:B}))
diff = sp.simplify(sp.expand(E0 - Lp6n))
print("(b) E_arc(kap=0, atoms->0) - 6*Lp(kap=0) =", diff)
print(f"({time.time()-t0:.0f}s)")
