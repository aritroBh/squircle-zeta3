"""thm:hexagon5 says the coefficients are polynomials in the edge data divided by L^5: list every
negative power that occurs in the 307 coefficients of the notes' table (read only)."""
import pickle, sympy as sp, collections
TAB = '<project>/research/08_theory/K1g_hexagon_order5_assembly/coef_acc_raw.pkl'
coef = pickle.load(open(TAB, 'rb'))
neg = collections.Counter(); syms = set()
for k, e in coef.items():
    syms |= e.free_symbols
    for p in e.atoms(sp.Pow):
        if p.exp.is_negative or (p.exp.is_Rational and not p.exp.is_integer):
            neg[(str(p.base), str(p.exp))] += 1
print('free symbols:', sorted(map(str, syms)))
print('negative or fractional powers (base, exponent): count of coefficients containing them')
for k, v in sorted(neg.items()):
    print('  ', k, v)
print('x-powers of the monomials:', sorted(collections.Counter(k[1] for k in coef).items()))
