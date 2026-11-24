import sys
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
d, m = T4.d, T4.m
for (dl,s) in [(7,1),(7,2),(8,1),(8,2),(8,3)]:
    f = T4.forced(dl,s)
    num, den = sp.fraction(sp.cancel(f))
    print(f"F({dl},{s}) = ( {sp.expand(num)} ) / ( {sp.factor(den)} )")
    print()
