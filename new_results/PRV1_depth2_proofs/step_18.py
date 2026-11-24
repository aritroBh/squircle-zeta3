
import sympy as s
from itertools import product, combinations_with_replacement, permutations

t = s.symbols('t', positive=True)
letters = [t, 1-t, 1+t, s.Integer(2)]

def logs(f):
    out = [0]*4
    num, den = s.fraction(s.factor(f))
    for a, sgn in [(num, 1), (den, -1)]:
        const, facs = s.factor_list(a, t)
        const = abs(const)
        out[3] += sgn*(s.factorint(s.numer(const)).get(2, 0)-s.factorint(s.denom(const)).get(2, 0))
        for fac, p in facs:
            j = next((j for j, z in enumerate(letters[:3]) if s.simplify(fac/z).is_number), None)
            if j is None: raise ValueError((f, fac))
            out[j] += sgn*p
    return s.Matrix(out)

def tens(*v):
    return s.Matrix([s.prod(v[j][i] for j, i in enumerate(q)) for q in product(range(4), repeat=len(v))])

def sy3(z): return -tens(logs(1-z), logs(z), logs(z))
def prod21(z, l):
    a, b, c = logs(1-z), logs(z), logs(l)
    return -tens(a,b,c)-tens(a,c,b)-tens(c,a,b)
def prod111(*l):
    v = list(map(logs, l))
    return sum((tens(*(v[j] for j in p)) for p in permutations(range(3))), s.zeros(64,1))

R,A,B = map(logs, letters[:3])
F_dec = -2*tens(R,A,A) + 2*tens(B,R,A) + tens(B,R,R)

# "good" z: at t=1/phi all values lie in {1/phi, -1/phi, 1/phi^2} and stay in (-1,1)
zlist = [t, -t, 1-t, t**2, 1-t**2, 1/(1+t), t/(1+t)]
products = [(z,l) for z in zlist for l in letters]
pure = list(combinations_with_replacement(letters,3))
cols = [sy3(z) for z in zlist]+[prod21(z,l) for z,l in products]+[prod111(*l) for l in pure]
labels = [s.polylog(3,z) for z in zlist]+[s.polylog(2,z)*s.log(l) for z,l in products]+[s.prod(s.log(v) for v in l) for l in pure]
M = s.Matrix.hstack(*cols)
print('symbol rank =', M.rank(), 'augmented rank =', M.row_join(F_dec).rank())
if M.rank() == M.row_join(F_dec).rank():
    sol, params = M.gauss_jordan_solve(F_dec)
    sol0 = sol.subs({p:0 for p in params})
    for n, v in zip(labels, sol0):
        if v: print(v, n)
