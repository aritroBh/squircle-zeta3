import pickle, sympy as sp
from sympy import Add, Function, sqrt, simplify, radsimp, expand
l2,l3,lu = sp.symbols('l2 l3 lu')
Li2A = Function('Li2A'); Li3A = Function('Li3A')
R = pickle.load(open('R256_sym.pkl','rb'))
# simplify each polylog argument
def simp_atoms(e):
    args = sorted(e.atoms(Li2A) | e.atoms(Li3A), key=str)
    mapping = {}
    for a in args:
        f = a.func; z = a.args[0]
        z2 = radsimp(simplify(z))
        mapping[a] = f(z2)
    return e.subs(mapping, simultaneous=True), mapping
R2, mapping = simp_atoms(R)
print("atom argument simplifications:")
for a,b in mapping.items():
    print("  ", a, " => ", b.args[0])
terms = Add.make_args(R2)
print("\nR(2,5,6) =")
for t in sorted(terms, key=str):
    print("   ", t)
