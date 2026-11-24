import sympy as sp
import mpmath as mp
c = sp.sqrt(3)/2
for el in [c, 1-c, 1+c]:
    print("el:", el, " has(c):", el.has(c), " subs:", el.subs({c: sp.Symbol('cvar')}))
# roots of 4x^2-3
roots = mp.polyroots([mp.mpf(4), mp.mpf(0), mp.mpf(-3)])
print("roots:", roots)
