
import sympy as s
from itertools import product, combinations_with_replacement, permutations

# Test: does G_c(t) = Li_21(-t, -(1-t)/(c t)) have a single-polylog antiderivative for c=3?
# Curve: x=-t, p=xy=(1-t)/c.  G_c' = Li11/t + Li2(-t)/(t(1-(1-c)t)) + c Li2((1-t)/c)/((1-t)(1-(1-c)t))
# For c=3: letters t, 1-t, 1+2t, 3. Compute symbol of G_3 and test membership in single-polylog span.
c_val = 3
t = s.symbols('t', positive=True)
letters = [t, 1-t, 1+2*t, s.Integer(3)]

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
    a, b, cc = logs(1-z), logs(z), logs(l)
    return -tens(a,b,cc)-tens(a,cc,b)-tens(cc,a,b)
def prod111(*l):
    v = list(map(logs,l))
    return sum((tens(*(v[j] for j in p)) for p in permutations(range(3))), s.zeros(64,1))

# Symbol of G_3 from dG = Li11 dlog(x) + [Li2(x)-Li2(p)/y] * dy/(1-y) structure.
# Simpler: dG = Li11 dlog t + Li2(-t) dlog(t(1+2t))/... let's use the derivative formula:
# G' = Li11/t + Li2(-t)/(t(1+2t)) + 3 Li2((1-t)/3)/((1-t)(1+2t))
# dG = Li11 dlog t + Li2(-t)[dlog t - 2 dlog(1+2t)] ... 1/(t(1+2t)) = 1/t - 2/(1+2t)
#     + Li2((1-t)/3)*3/((1-t)(1+2t));  3/((1-t)(1+2t)) = 1/(1-t) + 2/(1+2t)
# S(Li11): Li11 = log(1+t)log(1-y')... we need Li11(-t, -(1-t)/(3t)).
# z = -y(1-x)/(1-y): y=-(1-t)/(3t), 1-y=(1+2t)/(3t), 1-x=1+t => z = (1-t)(1+t)/(1+2t)?? 
#   = (1-t)(1+t)/(1+2t) = (1-t^2)/(1+2t)
# Li11 = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(z)
# -y/(1-y) = (1-t)/(1+2t)
R,A,B2 = map(logs, letters[:3])   # R=t, A=1-t, B2=1+2t
# S(log(1+t)log((1+2t)/(3t))): 1+t is NOT a letter! 1-x=1+t not in {t,1-t,1+2t,3}
try:
    logs(1+t)
    print("1+t factors OK")
except ValueError:
    print("1+t does NOT factor into letters {t,1-t,1+2t,3} -> Li11 itself outside class")
# So for c=3 even the Li11 closed form needs letters beyond the curve's alphabet (1+t appears).
# Check c=2 alphabet: letters {t,1-t,1+t,2}: 1-x = 1+t OK. c=1: 1+t OK (letters t,1-t,1+t).
# Conclusion: c=3 fails at the Li11 level already.
