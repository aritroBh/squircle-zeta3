
import sympy as sp
import sym_G
G21, G12, syms = sym_G.build()
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
sq = sp.sqrt(S0**2+Sc**2)/S0   # sqrt(1+Sc^2/S0^2)
# Remove the common sqrt factor and denominator: write G = sq * Num / Den
for label,G in [('G21',G21),('G12',G12)]:
    Ge = sp.expand(G).subs(sp.sqrt(1+Sc**2/S0**2), sq)
    Ge = sp.expand(Ge)
    num, den = sp.fraction(sp.together(Ge))
    num = sp.expand(num); den=sp.expand(den)
    print("====",label,"====")
    print("DEN =", sp.factor(den))
    print("NUM factored =", sp.factor(num))
    print()
