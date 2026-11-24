step = r'''
# step_02.py -- COR6X referee: controls (Task 1a, 1b) + deep consistency checks
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')

print("### Task 1(b): S(H_{0,0,0,1}(x)) = S(Li_4(x))")
S = sym_H((0,0,0,1), x)
print(S)
target = {((1-x), x, x, x): -1}
print("  equals -[1-x][x][x][x]:", sym_sub(S, target) == {})

print("\n### Li_n symbols (outermost integration letter last)")
for n in (2,3,4):
    w = (0,)*(n-1) + (1,)
    print(f"  S(Li_{n}(x)) =", sym_G(w, x))

print("\n### integrability of all symbols computed here")
for name, S in [("Li4", sym_H((0,0,0,1), x)), ("Li2", sym_G((0,1), x)),
                ("G(1,-1;x)", sym_G((1,-1), x)), ("Li22", sym_G((0,1/x,0,1/(x*y)), 1))]:
    r = integrable(S)
    print(f"  {name}: integrable = {r == {}}")

print("\n### Task 1(a): delta22 kills Li_4 and products")
SLi2x, SLi2y = sym_G((0,1), x), sym_G((0,1), y)
SLi3x = sym_G((0,0,1), x); Slogx = sym_G((0,), x)
print("  delta22(Li_4(x)) =", delta22(S))
print("  delta22(Li_2(x)^2) =", delta22(sym_shuffle(SLi2x, SLi2x)))
print("  delta22(Li_2(x)Li_2(y)) =", delta22(sym_shuffle(SLi2x, SLi2y)))
print("  delta22(log x * Li_3(x)) =", delta22(sym_shuffle(Slogx, SLi3x)))
print("  delta22(Li_1^4-type) =", delta22(sym_shuffle(sym_shuffle(SLi2x, SLi2x), sym_shuffle(SLi2x, SLi2x))))

print("\n### five-term relation for Li_2 modulo products (B_2)")
z1 = (1-x)/(1-x*y); z2 = (1-y)/(1-x*y)
terms = [sym_G((0,1), x), sym_G((0,1), y), sym_G((0,1), z1), sym_G((0,1), z2), sym_G((0,1), 1-x*y)]
tot = {}
for S in terms:
    tot = sym_add(tot, delta11(S))
print("  five-term (all +): ", tot)
if tot != {}:
    tot2 = {}
    for S in terms[:2] + [sym_scale(t,-1) for t in terms[2:]]:
        tot2 = sym_add(tot2, delta11(S))
    print("  variant (++---):", tot2)

print("\n### shuffle checks at total weight 2,3,4, letters {0,+-1} (eps-regularization test)")
import itertools
def check_shuffle(u, v):
    Su, Sv = sym_G(u, x), sym_G(v, x)
    lhs = sym_shuffle(Su, Sv)
    rhs = {}
    for w in shuffle_words(u, v):
        rhs = sym_add(rhs, sym_G(w, x))
    return sym_sub(lhs, rhs) == {}
words = [(0,), (1,), (-1,), (0,1), (0,-1), (1,-1), (-1,1), (0,0,1), (0,1,-1), (1,0,-1), (0,0,0,1)]
pairs = [((0,),(0,)), ((1,),(-1,)), ((0,),(1,)), ((0,1),(0,1)), ((0,),(0,1)), ((1,-1),(0,1)),
         ((0,1),(0,-1)), ((0,0,1),(0,)), ((0,0,1),(1,)), ((1,),(0,0,1)), ((0,),(0,0,0,1)),
         ((0,1),(1,-1)), ((1,), (1,0,-1))]
ok = True
for u, v in pairs:
    good = check_shuffle(u, v)
    ok = ok and good
    print(f"  shuffle {u} x {v}: {'OK' if good else 'FAIL'}")
print("  ALL SHUFFLE CHECKS:", "PASS" if ok else "FAIL")

print("\n### depth-two positive control: Li_{2,2}(x,y), and Li_{3,1}(x,y)")
S22 = sym_G((0, 1/x, 0, 1/(x*y)), 1)
print("  delta22(Li_{2,2}) nonzero:", delta22(S22) != {}, f"({len(delta22(S22))} terms)")
S31 = sym_scale(sym_G((0, 0, 1/x, 1/(x*y)), 1), -1)  # Li_{3,1} = -G(0,0,1/x,1/(xy);1)
print("  delta22(Li_{3,1}) nonzero:", delta22(S31) != {}, f"({len(delta22(S31))} terms)")
print("  Li22 integrable:", integrable(S22) == {}, " Li31 integrable:", integrable(S31) == {})
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_02.py","w").write(step)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_02.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout)
print("STDERR:", (r.stderr or "")[-1500:])
