
import sys
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)
W12, W1, W2 = ((0,1),(0,2)), ((0,1),(1,2)), ((0,2),(1,2))

# DGR three functions (planner note): Li22(-1,c), Li22(1/2, 2c/(1+c)), Li22(1/2, 2c/(c-1))
# letters:  Li22(x,y): X=[x],Y=[y],P=[1-x],Q=[1-y],R=[1-xy]
dgr1 = delta22(li22_sym(zero, cv, zero, mv, pv), 3)                       # Li22(-1, c)
dgr2 = delta22(li22_sym(zero, lsub(cv,pv), zero, lsub(mv,pv), lneg(pv)), 3)  # Li22(1/2, 2c/(1+c))
dgr3 = delta22(li22_sym(zero, lsub(cv,mv), zero, lsub(pv,mv), lneg(mv)), 3)  # Li22(1/2, 2c/(c-1))
def show(name, d):
    nm = {W12:'e1^e2', W1:'e1^e3', W2:'e2^e3'}
    print(f"{name}: " + ", ".join(f"{nm[k]}: {{{','.join(f'c^{e}:{v}' for e,v in sorted(p.items()))}}}" for k,p in sorted(d.items())))
show("Li22(-1,c)          ", dgr1)
show("Li22(1/2,2c/(1+c))  ", dgr2)
show("Li22(1/2,2c/(c-1))  ", dgr3)
target = {W1: {0: Fr(1), -1: Fr(-1)}, W2: {0: Fr(-1), -1: Fr(-1)}}
# solve target = b1*dgr1 + b2*dgr2 + b3*dgr3 with Laurent b's (unknowns per exponent -2..2)
import sympy as sp
unk = []
for j in range(3):
    for e in range(-2,3):
        unk.append((j,e))
rows = {}
for j, d in enumerate([dgr1,dgr2,dgr3]):
    for comp, poly in d.items():
        for e0, v0 in poly.items():
            for (jj,e) in unk:
                if jj != j: continue
                rows.setdefault((comp, e0+e), {})[(j,e)] = v0
tv = {}
for comp, poly in target.items():
    for e, v in poly.items():
        tv[(comp,e)] = v
A = sp.Matrix([[rows.get(key,{}).get(u, Fr(0)) for u in unk] for key in sorted(set(rows)|set(tv))])
b = sp.Matrix([tv.get(key, Fr(0)) for key in sorted(set(rows)|set(tv))])
sol = sp.linsolve((A,b), sp.symbols('x0:%d'%len(unk)))
print("DGR solve:", sol)
