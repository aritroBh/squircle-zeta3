"""M4e: factor-7 test for N5. All 4-subsets of the 8 Li32 candidates;
solve the cobracket correction, report denominator factors.
"""
import json, itertools, sys
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import flint
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce as cred, target as ctarget, depth2
from w4_right_coproduct import ZERO, ONE, Q, R, A, B

HERE = Path(__file__).resolve().parent
def log(s): print(s, flush=True)
def main():
    C.ban_float()
    lhs = ctarget()
    log(f'target D = {lhs}')
    cand = [(B,ONE),(R/2,ONE),(-B,ONE),(B,-ONE),(R/2,-ONE),(ONE-R/2,ONE),(2-R,ONE),(-B*B,ONE)]
    reds = [depth2(x,y) for x,y in cand]
    atoms = [('z3',),(('Li3',B),),(('Li3',-B),),(('Li3',R/2),)]
    rhs = [Q(F(str(lhs.d.get(at,ZERO).p)), F(str(lhs.d.get(at,ZERO).q))) for at in atoms]
    M = [[Q(F(str(red.d.get(at,ZERO).p)), F(str(red.d.get(at,ZERO).q))) for at in atoms] for red in reds]
    def qdet(rows4):
        # 4x4 determinant over Q(sqrt2) via sympy
        import sympy as sp
        m = sp.Matrix([[c.p + c.q*sp.sqrt(2) for c in row] for row in rows4])
        return sp.nsimplify(m.det())
    results = []
    for idxs in itertools.combinations(range(8), 4):
        rows4 = [M[i] for i in idxs]
        # solve rows4^T x = rhs  (x_i multiplies cobracket of atom i)
        import sympy as sp
        mm = sp.Matrix([[rows4[r][c].p + rows4[r][c].q*sp.sqrt(2) for c in range(4)] for r in range(4)])
        if mm.det() == 0:
            continue
        rv = sp.Matrix([c.p + c.q*sp.sqrt(2) for c in rhs])
        sol = mm.inv()*rv  # sol[c] = coefficient for atom index idxs[c]... careful transpose
        # we need sum_i x_i M[idxs[i]][atom] = rhs[atom]: x = (M4)^{-1} rhs with M4 rows=atoms? redo:
        M4 = sp.Matrix([[M[i][a].p + M[i][a].q*sp.sqrt(2) for i in idxs] for a in range(4)])  # rows=atoms, cols=sel
        if M4.det() == 0: continue
        x = M4.inv()*sp.Matrix([c.p + c.q*sp.sqrt(2) for c in rhs])
        # denominators
        dens = []
        for v in x:
            v = sp.nsimplify(v)
            rp = sp.fraction(sp.nsimplify(v).as_numer_denom()[1]) if False else None
        # get denominator over Z: write v = u + w sqrt2, lcm of denoms
        alld = 1
        for v in x:
            v = sp.expand(v)
            p = sp.Rational(v.coeff(sp.sqrt(2), 0))
            q = sp.Rational(v.coeff(sp.sqrt(2), 1))
            alld = lcm(alld, p.denominator, q.denominator)
        fac = sp.factorint(alld)
        results.append((alld, idxs, fac))
    results.sort()
    log('denominator | factors | atom 4-subset (smallest 12):')
    for alld, idxs, fac in results[:12]:
        names = [f'({cand[i][0]},{cand[i][1]})' for i in idxs]
        log(f'  den={alld} fac={fac} {names}')
    sevenfree = [r for r in results if 7 not in r[2]]
    log(f'7-free subsets: {len(sevenfree)} of {len(results)}')

if __name__=='__main__':
    main()
