# step_23: generic exact vecfun for Q(sqrt2),Q(sqrt3),Q(sqrt5) + 12-gon & decagon CR sets
import sympy as sp
from fractions import Fraction as Fr
import mpmath as mp

class QFVec:
    def __init__(self, D, letters):
        """letters: list of (name, elem); assumes every z in sight is +-prod letters^a.
        Solve via norm (prime content) + log|z|, then VERIFY exactly."""
        self.D = D; self.sq = sp.sqrt(D)
        self.letters = letters
        self.norms = [sp.simplify(e*e.subs(self.sq,-self.sq)) for _,e in letters]
    def vec(self, z):
        z = sp.simplify(sp.expand(z))
        N = sp.simplify(z*z.subs(self.sq,-self.sq))
        N = sp.Rational(N)
        # solve for exponents: |N| = prod |Ni|^ai (prime factorization), then logs for the rest
        aN = abs(N)
        # first: prime-content solve for letters with nontrivial norm
        import itertools
        # build exponent vector by linear algebra over prime exponents of |N|
        primes = set()
        facN = sp.factorint(aN.numerator); facNd = sp.factorint(aN.denominator)
        for p in list(facN)+list(facNd): primes.add(p)
        lfacs = []
        for ni in self.norms:
            ani = abs(sp.Rational(ni))
            fn = sp.factorint(ani.numerator); fd = sp.factorint(ani.denominator)
            for p in list(fn)+list(fd): primes.add(p)
        primes = sorted(primes)
        def pvec(q):
            q = abs(sp.Rational(q)); fn = sp.factorint(q.numerator); fd = sp.factorint(q.denominator)
            return [fn.get(p,0)-fd.get(p,0) for p in primes]
        A = sp.Matrix.hstack(*[sp.Matrix(pvec(ni)) for ni in self.norms])
        b = sp.Matrix(pvec(N))
        sol = A.gauss_jordan_solve(b)[0]
        exps = [sp.Rational(sol[i,0].subs({s_:0 for s_ in sol.free_symbols})) for i in range(len(self.letters))]
        # residual unit part: z / prod e_i^exps should be a unit +-1... if not, fix by log
        prod = sp.Integer(1)
        for (_,e),a in zip(self.letters, exps):
            prod *= e**a
        rem = sp.simplify(sp.expand(z) - sp.expand(prod))
        rem2 = sp.simplify(sp.expand(z) + sp.expand(prod))
        if rem == 0 or rem2 == 0:
            return {self.letters[i][0]: Fr(exps[i].numerator, exps[i].denominator) for i in range(len(exps)) if exps[i]}
        # fix via log
        mp.mp.dps = 60
        ratio = sp.simplify(sp.expand(z)/sp.expand(prod))
        lr = mp.log(abs(mp.mpf(str(sp.N(ratio,50)))))
        # try adjusting the last letters (units) by brute force small integers
        import itertools
        for delta in itertools.product(range(-6,7), repeat=len(exps)):
            if sum(d!=0 for d in delta)==0: continue
            prod2 = prod
            for (_,e),d in zip(self.letters, delta):
                prod2 *= e**d
            if sp.simplify(sp.expand(z)-sp.expand(prod2))==0 or sp.simplify(sp.expand(z)+sp.expand(prod2))==0:
                ex = [exps[i]+delta[i] for i in range(len(exps))]
                return {self.letters[i][0]: Fr(ex[i].numerator, ex[i].denominator) for i in range(len(ex)) if ex[i]}
        raise ValueError(f"vec failed for {z} (N={N}, exps={exps})")

g2 = 1+sp.sqrt(3); u3 = 2+sp.sqrt(3)
V3 = QFVec(3, [("g2",g2),("g3",sp.sqrt(3)),("u",u3)])
# verify on table entries
T = pickle.load(open("field_tables.pkl","rb")) if False else None
tests3 = [2-sp.sqrt(3),(sp.sqrt(3)-1)/2,sp.sqrt(3)/3,7-4*sp.sqrt(3),(2-sp.sqrt(3))/2,(2*sp.sqrt(3)-3)/3,(3-sp.sqrt(3))/6,sp.Rational(1,3),
          sp.sqrt(3)-1,(3-sp.sqrt(3))/2,(3-sp.sqrt(3))/3,-6+4*sp.sqrt(3),sp.sqrt(3)/2,(6-2*sp.sqrt(3))/3,(3+sp.sqrt(3))/6,sp.Rational(2,3),
          (3+sp.sqrt(3))/2,4+2*sp.sqrt(3),-(3+2*sp.sqrt(3)),-(1+sp.sqrt(3))/2]
for z in tests3:
    print(z, "->", V3.vec(z))
