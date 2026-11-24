import sympy as sp, time, pickle, os
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
import sys; sys.path.insert(0, work)
d, m = sp.symbols('d m')
i = sp.I
# Recompute total(5) quickly (2s) and analyze
exec(open(os.path.join(work,'prove_S5.py').read().split('t0 = time.time()')[0].split('# ---- independent')[0].replace('import sympy as sp, time, pickle, os','').replace('_HERE = os.path.dirname(os.path.abspath(__file__))','_HERE="x"')))
T5 = total(5)
R5 = sp.together(T5)
num5, den5 = sp.fraction(R5)
print("denominator of total(5) factored:")
print(sp.factor(den5))
print("numerator degree:", sp.Poly(sp.expand(num5), d, m).total_degree() if hasattr(sp.Poly(sp.expand(num5),d,m),'total_degree') else sp.degree(sp.expand(num5)))
print(f"({round(time.time()-t0,1)}s)")
