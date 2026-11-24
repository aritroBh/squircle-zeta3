
import sys
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)
def lsub(u,v): return ladd(u, lneg(v))
W12, W1, W2 = ((0,1),(0,2)), ((0,1),(1,2)), ((0,2),(1,2))
dgr1 = delta22(li22_sym(zero, cv, zero, mv, pv), 3)                       # Li22(-1, c)
dgr2 = delta22(li22_sym(zero, lsub(cv,pv), zero, lsub(mv,pv), lneg(pv)), 3)  # Li22(1/2, 2c/(1+c))
dgr3 = delta22(li22_sym(zero, lsub(cv,mv), zero, lsub(pv,mv), lneg(mv)), 3)  # Li22(1/2, 2c/(c-1))
def show(name, d):
    nm = {W12:'e1^e2', W1:'e1^e3', W2:'e2^e3'}
    print(f"{name}: " + ", ".join(f"{nm[k]}: {{{','.join(f'c^{e}:{v}' for e,v in sorted(p.items()))}}}" for k,p in sorted(d.items())))
show("Li22(-1,c)          ", dgr1)
show("Li22(1/2,2c/(1+c))  ", dgr2)
show("Li22(1/2,2c/(c-1))  ", dgr3)
# also my two:
show("Li22(1/2,1-c)       ", delta22(li22_sym(zero, mv, zero, cv, pv), 3))
show("Li31(1/2,1+c)       ", delta22(li31_sym(zero, pv, zero, cv, mv), 3))
show("Li22(1/2,1+c)       ", delta22(li22_sym(zero, pv, zero, cv, mv), 3))
show("Li31(1/2,1-c)       ", delta22(li31_sym(zero, mv, zero, cv, pv), 3))
