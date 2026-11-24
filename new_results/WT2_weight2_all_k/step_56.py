import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from fastred import v
S2 = sp.sqrt(2); a = 1+S2; b = S2-1; I = sp.I
eps = sp.Symbol('eps')

g20   = (1+S2)/4 - S2/(2*(v**2+1))
g11_3 = -S2*(v**2+2*v-1)/(4*(v**2+1))
R10   = (1-S2/2)/v - 2/(v+a)
Rw2   = 1/(2*v) + 1/(2*(v-a)) - 1/(v+a) - 1/(2*(v+b))
DL1 = 1/(v-b) + 1/(v+a) - 1/v
DL2 = 1/(v-a) + 1/(v+b) - 1/v

def F(R):   # v -> 1/v fold
    return sp.together(R.subs(v, 1/v)/v**2)

Gm20 = eps**2/(1-eps) * g20
Gm11 = -eps**2 * g20 + eps**3/(1-eps) * g11_3

Wcas = -2*Gm20*DL1 - Gm11*DL2 - F(sp.together(Gm11*DL1))
Wcas = sp.together(Wcas)

# partial fraction over Q(sqrt2, I): check no poles at +-i, and value
pf = sp.apart(Wcas, v, extension=[S2, I])
terms = sp.Add.make_args(sp.expand(pf))
ipolet = [t for t in terms if sp.denom(sp.together(t)).has(I)]
real = sp.together(sum(t for t in terms if not sp.denom(sp.together(t)).has(I)))
print("poles at +-i in Wcas:", ipolet if ipolet else "NONE (cancel)")
target = eps**2/(1-eps)*Rw2
print("SimplePart(Wcas) == eps^2/(1-eps) Rw2 :", sp.simplify(sp.radsimp(sp.together(real - target)))==0)
print()
print("=> S_10 = eps*R10 + eps^2/(1-eps)*Rw2  [conditional on Gm20/Gm11 stabilization]")
print("=> [eps^k] S_10 = Rw2 for all k>=2")
