"""Numeric faithfulness check: the G-expression built (by M5a's gexpr5.g_expression) from a fed kernel vector equals the
direct quadrature of  int_b^a sum_{ij,r} e_{ij,r} log^i t log^j(1-t) dv/(v-r).  Usage: numeric_check.py E55"""
import sys, json, time
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT/'orig'))
import gexpr5 as G                       # ORIGINAL M5a g_expression builder
import geval_lib as E                    # ORIGINAL M5a numeric G evaluator (lower precision)
import mpmath as mp
mp.mp.dps = 60
name = sys.argv[1]
om = json.load(open(ROOT/'cases'/name/'omega5.json'))
parts = []
for ij in [(3,2),(4,1),(5,0)]:
    parts.append(G.g_expression(ij[0], ij[1], om[f"{ij[0]},{ij[1]}"]))
gex = G.combine(*parts)
tsv = ROOT/'cases'/name/'gpl.tsv'
G.save_tsv(gex, str(tsv), "# gexpr5.g_expression of fed kernel")
print(name, 'G-terms', len(gex), flush=True)
t0 = time.time()
val, nreg = E.evaluate(str(tsv))
print(f'G-expression evaluated in {time.time()-t0:.0f}s, endpoint-regularized words {nreg}, flint rad {val.real.rad().str(5)}, imag {val.imag.mid().str(5)}', flush=True)
got = mp.mpf(val.real.mid().str(55, radius=False, more=True))
# direct quadrature of the kernel integral (independent of every G-function)
r2 = mp.sqrt(2); a = 1+r2; b = r2-1
pole = {'0': mp.mpf(0), 'a': a, '-a': -a, '-b': -b}
cells = []
for ij in [(3,2),(4,1),(5,0)]:
    for nm, (p, q) in om[f"{ij[0]},{ij[1]}"]:
        cells.append((ij[0], ij[1], pole[nm], mp.mpf(F(p).numerator)/F(p).denominator + r2*mp.mpf(F(q).numerator)/F(q).denominator))
def integrand(v):
    t = (v*v+2*v-1)/(4*v)
    lt, lm = mp.log(t), mp.log(1-t)
    return sum(c*lt**i*lm**j/(v-z) for i, j, z, c in cells)
quad = mp.quad(integrand, [b, mp.mpf(1)/2, 1, 2, a], maxdegree=10)
print('G-expression value =', mp.nstr(got, 45))
print('direct quadrature  =', mp.nstr(quad, 45))
print('|difference|       =', mp.nstr(abs(got-quad), 5), flush=True)
