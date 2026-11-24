step = r'''
# step_02.py -- COR6X referee: controls (Task 1a, 1b) + deep consistency checks
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp
x, y = sp.symbols('x y')

print("### Task 1(b): S(H_{0,0,0,1}(x)) = S(Li_4(x))")
SLi4f = sym_H((0,0,0,1), x)
target = {((1-x), x, x, x): -1}
print("  S(H_{0,0,0,1})(x) =", SLi4f, " equals -[1-x][x][x][x]:", sym_sub(SLi4f, target) == {})

print("### S(Li_n(x)) = -[1-x][x]...[x] (outermost integration letter last)")
for n in (2,3,4):
    SLi = sym_scale(sym_G((0,)*(n-1) + (1,), x), -1)
    print(f"  Li_{n}:", sym_sub(SLi, {(1-x,) + (x,)*(n-1): -1}) == {})

print("### integrability (differential 2-form criterion)")
for name, S, vs in [("Li4", sym_H((0,0,0,1), x), (x,)), ("Li2", sym_G((0,1), x), (x,)),
                    ("G(1,-1)", sym_G((1,-1), x), (x,)),
                    ("Li22", sym_G((0,1/x,0,1/(x*y)), 1), (x,y)),
                    ("Li31", sym_G((0,0,1/x,1/(x*y)), 1), (x,y))]:
    print(f"  {name}: integrable = {integrable(S, vs) == {}}")

print("### Task 1(a): delta22 kills Li_4 and all products (expanded letters)")
SLi2x, SLi2y = sym_G((0,1), x), sym_G((0,1), y)
SLi3x = sym_G((0,0,1), x); Slogx = sym_G((0,), x); Slog1mx = sym_G((1,), x)
l1 = sym_shuffle(Slogx, Slog1mx)
for name, S in [("Li_4(x)", SLi4f), ("Li_2(x)^2", sym_shuffle(SLi2x, SLi2x)),
                ("Li_2(x)Li_2(y)", sym_shuffle(SLi2x, SLi2y)),
                ("log x * Li_3(x)", sym_shuffle(Slogx, SLi3x)),
                ("logx log(1-x) Li_2(x)", sym_shuffle(l1, SLi2x)),
                ("log^2x log^2(1-x)", sym_shuffle(l1, l1))]:
    print(f"  delta22({name}) =", delta22(sym_expanded(S)))

print("### five-term relation for Li_2 (expanded): {x}+{y}+{z1}+{z2}+{1-xy} = 0")
z1 = (1-x)/(1-x*y); z2 = (1-y)/(1-x*y)
tot = {}
for z in (x, y, z1, z2, 1-x*y):
    tot = sym_add(tot, delta11(sym_expanded(sym_G((0,1), z))))
print("  =", tot)

print("### shuffle checks, letters {0,+-1} (eps-regularization test)")
def check_shuffle(u, v):
    lhs = sym_shuffle(sym_G(u, x), sym_G(v, x))
    rhs = {}
    for w in shuffle_words(u, v):
        rhs = sym_add(rhs, sym_G(w, x))
    return sym_sub(lhs, rhs) == {}
pairs = [((0,),(0,)), ((1,),(-1,)), ((0,),(1,)), ((0,1),(0,1)), ((0,),(0,1)), ((1,-1),(0,1)),
         ((0,1),(0,-1)), ((0,0,1),(0,)), ((0,0,1),(1,)), ((1,),(0,0,1)), ((0,),(0,0,0,1)),
         ((0,1),(1,-1)), ((1,),(1,0,-1)), ((-1,),(0,1,-1)), ((0,-1),(0,0,1))]
ok = all(check_shuffle(u, v) for u, v in pairs)
print("  ALL", len(pairs), "SHUFFLE CHECKS:", "PASS" if ok else "FAIL")

print("### depth-two positive controls (expanded)")
S22 = sym_G((0, 1/x, 0, 1/(x*y)), 1)
D22 = delta22(sym_expanded(S22))
print("  delta22(Li_{2,2}(x,y)):", "NONZERO" if D22 else "zero", f"({len(D22)} terms)")
for k, v in sorted(D22.items(), key=str): print("    ", k, v)
S31 = sym_scale(sym_G((0, 0, 1/x, 1/(x*y)), 1), -1)
D31 = delta22(sym_expanded(S31))
print("  delta22(Li_{3,1}(x,y)):", "NONZERO" if D31 else "zero", f"({len(D31)} terms)")
for k, v in sorted(D31.items(), key=str): print("    ", k, v)
# compare D22 with [Li_2(x)]^[Li_2(y)] = ((1-x)^x)^((1-y)^y):
w1 = ((1-x), x); w2 = ((1-y), y)
print("  equals ((1-x)^x)^((1-y)^y)?", D22 == {(w1, w2): 1} or D22 == {(w2, w1): -1})
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_02.py","w").write(step)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_02.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols")
print(r.stdout); print("STDERR:", (r.stderr or "")[-1500:])
open("<project>/research/08_theory/COR6X_independent_corner_symbols/progress.log","a").write("step_02 done: controls pass (Li_n, five-term, shuffles, delta22 controls)\n")
