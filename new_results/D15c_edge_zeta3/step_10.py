
# Verify total degree in (L0,L1) is <=3, so G21,G12 are the complete cross coefficients
degs = set()
for term in D3.as_ordered_terms():
    d0 = sp.degree(term, L0); d1 = sp.degree(term, L1)
    degs.add((d0,d1))
print("All (deg L0, deg L1) pairs present:", sorted(degs))
print("Max total degree:", max(a+b for a,b in degs))
print("Any term with degL0>=2 and degL1>=1 beyond (2,1)?", [(a,b) for a,b in degs if a>=2 and b>=1])
print("Any term with degL0>=1 and degL1>=2 beyond (1,2)?", [(a,b) for a,b in degs if a>=1 and b>=2])

# Simplify G21, G12
G21s = sp.factor(sp.simplify(G21))
G12s = sp.factor(sp.simplify(G12))
print("\nG21 simplified:")
print(G21s)
print("\nG12 simplified:")
print(G12s)
