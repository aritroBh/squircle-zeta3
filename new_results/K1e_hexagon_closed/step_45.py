
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB,rr = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB r')
lines=[]
for name in ['zeta3','Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','zeta2*ln(1-x)','zeta2','Li2(x)','ln^3(1-x)','ln^2(1-x)','ln(1-x)','1']:
    num,den = clean[name]
    n2 = sp.expand(num.subs({kap: rr-1, A1: 1-B}))
    # collect by log monomials
    p = sp.Poly(n2, lA0,lA1,lB)
    lines.append(f"===== coefficient of {name}  (denominator L2^(5/2); A0=1-B*r, A1=1-B, kap=r-1, x=1-r) =====")
    for exps,coef in sorted(p.terms(), reverse=True):
        mu = lA0**exps[0]*lA1**exps[1]*lB**exps[2]
        lines.append(f"  [{mu}] * ( {sp.factor(coef)} )")
    lines.append("")
txt="\n".join(lines)
open(f"{WORK}/formula_table.txt","w").write(txt)
print(txt[:3200])
print("... written to formula_table.txt, total chars:", len(txt))
