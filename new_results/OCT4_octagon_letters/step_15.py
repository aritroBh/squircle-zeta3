import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import log
import sympy as sp

sqrt2 = sp.sqrt(2)
# Normals n_i = (cos(i pi/4), sin(i pi/4)), i=1..4
normals = []
for i in range(1,5):
    th = sp.pi*i/4
    normals.append((sp.cos(th), sp.sin(th)))
print("normals:", normals)

# Vertices of K_1: solve pairs of active constraints <n_i,x>=±1 region.
# From analysis: top vertex (0, sqrt2-1), upper-right vertex (1-sqrt2/2, 1-sqrt2/2)
V0 = sp.Matrix([0, sqrt2-1])
V1 = sp.Matrix([1-sqrt2/2, 1-sqrt2/2])
print("V0 =", V0.T, " V1 =", V1.T)

# Support S(x) = sum |<n_i,x>|; check both vertices give S=1
def S(x):
    return sum(abs(sp.Matrix(n).dot(x)) for n in normals)
# symbolic: at V0
for name,V in [("V0",V0),("V1",V1)]:
    s = sp.simplify(sum(sp.Abs(sp.Matrix(n).dot(V)) for n in normals))
    print(name, "S =", sp.simplify(s))
log("step: octagon edge V0=(0,sqrt2-1), V1=(1-sqrt2/2,1-sqrt2/2); both on S=1")
