
import sympy as sp, sys, time
sys.path.insert(0, r"research/08_theory/D17_zeta5_polygons")
sys.path.insert(0, r"research/08_theory/D15c_edge_zeta3")
import sym_D5
import sym_G   # from D15c

# Build reduced D3
D, syms = sym_D5.build_D(3, reduced=True)
L0,L1 = syms['L0'], syms['L1']
D3 = sp.expand(D[3])
G21_red = sp.expand(D3.coeff(L0,2).coeff(L1,1))
G12_red = sp.expand(D3.coeff(L0,1).coeff(L1,2))
print("reduced G21 terms:", len(G21_red.as_ordered_terms()))

# sym_G G21,G12 (proved tau/gam-independent)
G21g, G12g, gsyms = sym_G.build()
# sym_G uses s0,s1,c0,c1,S0,Sc same names; D0 factor: sym_G D3 includes D0=sqrt(B)/S0
# Compare
d21 = sp.simplify(G21_red - G21g)
d12 = sp.simplify(G12_red - G12g)
print("G21_red - G21g =", d21)
print("G12_red - G12g =", d12)
