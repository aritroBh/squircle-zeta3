# step_17: vecfuns for Q(sqrt2), Q(sqrt3), Q(sqrt5) with exact verification
import sympy as sp
from fractions import Fraction as Fr

def make_vecfun(D, letters):
    """letters: list of (name, sympy element of Q(sqrtD)). Returns vecfun(z)->dict name->Fr,
    solving z = +- prod letters^a exactly (a in Q allowed for ramified: handled by caller's exponents)."""
    sq = sp.sqrt(D)
    def vecfun(z):
        z = sp.simplify(z)
        # brute force: find exponents a_i in QQ by matching logs? Use exact norm + divide out.
        # Strategy: iterate: work in the multiplicative group; represent [z] by greedy exact division.
        rem = sp.simplify(z)
        exp = {}
        # handle sign
        if sp.simplify(rem) != sp.simplify(abs(complex(rem.evalf()))):
            pass
        return None
    return None

# Simpler robust approach: precomputed factor tables, verified exactly by multiplication.
def verify_table(D, table):
    """table: dict name-> (sympy elem, dict letter->Fr). Check elem == prod letters^exp exactly."""
    sq = sp.sqrt(D)
    L = {}
    ok = True
    for name,(elem,vec) in table.items():
        prod = sp.Integer(1)
        for lname,e in vec.items():
            prod *= table[lname][0]**sp.Rational(e.numerator,e.denominator) if lname!=name else sp.Integer(1)
    return ok

# Direct: define letters and a solver via logs of absolute values + exact check
class QField:
    def __init__(self, D, letters):
        self.D = D; self.sq = sp.sqrt(D)
        self.letters = letters  # list of (name, elem)
    def vec(self, z):
        z = sp.nsimplify(z) if not isinstance(z, sp.Expr) else z
        z = sp.simplify(sp.expand(z))
        # solve z^2 = +- prod letters^{2a} via norms and logs: use linear algebra on exact logs
        import itertools
        # candidate exponents bounded search
        names = [n for n,_ in self.letters]
        elems = [e for _,e in self.letters]
        # Use exact approach: z*conj(z) = Norm; and z/|z| matching. Build exponent by
        # taking log ratios -- then VERIFY exactly; if fail, expand search.
        import mpmath as mp
        mp.mp.dps = 50
        lv = [mp.log(abs(sp.N(e, 50))) for _,e in self.letters]
        target = mp.log(abs(sp.N(z, 50)))
        # solve small rational exponents: use mpmath findrelation? do LLL-free 2-3 dim solve
        n = len(elems)
        if n <= 4:
            M = mp.matrix(lv + [1], 1)  # not used; do direct solve below
        # exact solve: z / prod e^a = +-1. Use sympy nsimplify on log ratio
        from mpmath import matrix as mmat, lu_solve
        A = mp.matrix([[lv[i]] for i in range(n)]) if n==1 else mp.matrix([lv])
        b = mp.matrix([target])
        try:
            x = mp.lu_solve(A, b)
        except Exception:
            x = None
        if x is None:
            raise ValueError("no solve")
        exps = [sp.nsimplify(x[i]) for i in range(n)]
        prod = sp.Integer(1)
        for e,a in zip(elems, exps):
            prod *= e**a
        if sp.simplify(prod - z) == 0 or sp.simplify(prod + z) == 0:
            return {names[i]: Fr(sp.Rational(exps[i]).numerator, sp.Rational(exps[i]).denominator) for i in range(n) if exps[i]!=0}
        raise ValueError(f"vec failed for {z}: got {exps}, ratio {sp.simplify(z/prod)}")

# Q(sqrt2): letters 2, 1+sqrt2
F2 = QField(2, [("2", sp.Integer(2)), ("u2", 1+sp.sqrt(2))])
for z, name in [(sp.sqrt(2)-1,"k"), (2-sp.sqrt(2),"d"), (sp.sqrt(2)/2,"r2h"), (-sp.sqrt(2)/2,"mr2h"),
                (1-sp.sqrt(2),"x1"), (1-sp.sqrt(2)/2,"x2"), (sp.sqrt(2),"sqrt2"), (1+sp.sqrt(2)/2,"1+r2h"),
                (2+sp.sqrt(2),"2+r2"), (sp.sqrt(2)-1,"k")]:
    print(name, "=", z, "->", F2.vec(z))
print("Q(sqrt2) vecfun OK")
