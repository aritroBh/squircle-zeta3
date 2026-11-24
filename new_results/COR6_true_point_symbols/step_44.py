from fractions import Fraction as Fr
import importlib, cor6_nf
importlib.reload(cor6_nf)
from cor6_nf import NF, primes_above

def show_target(label, minpoly, c0):
    F = NF(minpoly)
    c = c0
    one = F.one()
    two = F.from_int(2)
    Nc = F.norm(c); N1mc = F.norm(F.sub(one, c)); N1pc = F.norm(F.add(one, c))
    # S0 = primes in denominators/numerators of norms
    import sympy as sp
    S0 = set([2])
    for N in (Nc, N1mc, N1pc):
        S0 |= set(sp.factorint(N.numerator).keys()) | set(sp.factorint(N.denominator).keys())
    info = {}
    for p in sorted(S0):
        try:
            cnt, degs = primes_above(F.m, p)
        except Exception as e:
            cnt, degs = ("ERR", str(e))
        info[p] = (cnt, degs)
    print(f"{label}: deg {F.d}, N(c)={Nc}, N(1-c)={N1mc}, N(1+c)={N1pc}, S0={sorted(S0)}")
    for p in sorted(S0):
        print(f"    p={p}: #primes={info[p]}")

show_target("c=sqrt2/2 (square,n=4)", [-2, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sqrt3/2 (triangle,n=3)", [-3, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sin(pi/10) (n=10)", [-5, 0, 1], (Fr(-1, 4), Fr(1, 4)))
show_target("c=sin(pi/5) (n=5)", [5, 0, -5, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sin(pi/8)=cos(3pi/8) (n=8)", [2, 0, -4, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=sin(pi/12) (n=12)", [1, 0, -4, 0, 1], (Fr(0), Fr(1, 2)))
show_target("c=cos(pi/12) (alpha=pi/6)", [1, 0, -4, 0, 1], (Fr(0), Fr(1, 2)))  # same field, c0 = cos(pi/12) as poly in theta=2sin(pi/12)?
