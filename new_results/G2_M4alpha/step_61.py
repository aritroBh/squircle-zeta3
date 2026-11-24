import pickle, sympy as sp
cwd = "<project>/research/08_theory/G2_M4alpha"
n = sp.symbols('n')
terms = pickle.load(open(cwd + '/mu4_formula.pkl', 'rb'))
facs = set()
dens = set()
for mon, cf in terms:
    e2, e3, e4 = mon[8], mon[9], mon[10]
    if e2 or e3 or e4:
        continue  # only rational part
    num, den = sp.fraction(sp.together(sp.sympify(cf)))
    for f, mult in sp.factor_list(den)[1]:
        facs.add(str(sp.factor(f)))
print("denominator factors in the rational part:", sorted(facs))
# also for z2 part
facs2 = set()
for mon, cf in terms:
    if mon[8] == 1:
        num, den = sp.fraction(sp.together(sp.sympify(cf)))
        for f, mult in sp.factor_list(den)[1]:
            facs2.add(str(sp.factor(f)))
print("denominator factors in the z2 part:", sorted(facs2))
