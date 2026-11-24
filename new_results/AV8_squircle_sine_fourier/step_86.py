from mpmath import mp, mpf, pi, log, zeta, pslq, nstr
mp.dps = 50
import json
res = json.load(open('scratchpad/bnder60.json'))
ln2 = log(2)
def bp(n):
    n=int(n)
    if n==1: return (ln2-1)/2
    return mpf((-1)**((n-1)//2))/((n**2-1)*(2*((n-1)//4)+1))
# measured b_n''
bpp = {n: mpf(res[f"{n},2"]) for n in [1,3,5,7,9,11,13]}
S13 = sum(bp(n)*bpp[n] for n in [1,3,5,7,9,11,13])
print("partial Sigma (n<=13):", nstr(S13,40))
# tail estimate: b_n'' asymp ~ ? use C_n law + B_n=b_n' + A_n unknown; rough magnitude of tail:
# |b_n' b_n''| ~ (2/n^3)*(~1/(2n^2)) = 1/n^5, tail > 13 ~ 1/(4*13^4) ~ 8.8e-6
import sympy as sp
p = sp.symbols('p', positive=True)
pip = 2*sp.pi/(p*sp.sin(sp.pi/p))
N = 4/(p*pip)*sp.beta(3/p, 1-1/p)
N3 = sp.diff(N, p, 3).subs(p, 2)
N3 = N3.replace(lambda e: isinstance(e, sp.log) and e.args[0].is_Integer,
                lambda e: sum(m*sp.log(q) for q,m in sp.factorint(e.args[0]).items()))
N3 = sp.simplify(N3.subs({sp.polygamma(1, sp.Rational(3,2)): sp.pi**2/2-4,
                          sp.polygamma(2, sp.Rational(3,2)): -2*sp.zeta(3)+12}))
print("N'''(2) =", N3)
print("float:", sp.N(N3,30))
