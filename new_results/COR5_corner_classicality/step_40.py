import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
import cor5_sym
from cor5_sym import W4Quotient, li_sym, rref
sols = json.load(open("sunit23_solutions.json"))
def vec23(fac):
    s, a, b = fac
    return (Fr(a), Fr(b))
Q = W4Quotient(2)
li4 = [(Fr(xs), Q.coords(li_sym(4, vec23(fx), vec23(fy)))) for xs, fx, fy in sols]
vF = (Fr(5,2), Fr(1,2), Fr(0))
M = [list(v) for _, v in li4] + [list(vF)]
out = rref(M)
print("module rref rank:", len(out[0]), "pivots:", out[1])
import inspect
print(inspect.getsource(cor5_sym.rref))
