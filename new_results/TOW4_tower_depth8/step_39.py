import sys
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import sympy as sp
import tow2lib as T2
i = sp.I
# compute tau1(8), tau2(8) from definitions, simplify into Q(i)
for n in [8]:
    t1 = T2.tau1(n); t2 = T2.tau2(n)
    # reduce sqrt2 algebraically: express as A + B*sqrt2 then assert B==0
    t1e = sp.expand(t1).rewrite(sp.exp) if False else t1
    t1r = sp.simplify(sp.expand(t1))
    t2r = sp.simplify(sp.expand(t2))
    print("tau1(8) raw simplified:", t1r)
    print("tau2(8) raw simplified:", t2r)
    # algebraic reduction: separate sqrt2 part
    t1_alg = sp.algebraic_field if False else None
    # Use eval at sqrt2 symbol: collect
    s2 = sp.sqrt(2)
    t1c = sp.collect(sp.expand(t1), s2)
    t2c = sp.collect(sp.expand(t2), s2)
    print("tau1(8) collected:", t1c)
    print("tau2(8) collected:", t2c)
