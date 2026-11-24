import json, mpmath, itertools
mpmath.mp.dps = 405
from mpmath import mpf, log, pi, zeta
from flint import fmpz_mat, fmpz
lam = log(2)
# basis: monomials lam^a pi^(2b) prod zeta(odd)^mult, weight w
def monoms(w):
    out = []
    def rec(rem, factors, start):
        if rem == 0:
            out.append(list(factors)); return
        for g in range(start, rem+1):
            if g in (1,2) or g%2==1:
                rec(rem-g, factors+[g], g)
    rec(w, [], 1)
    return out
def mval(m):
    v = mpf(1)
    for g in m:
        v *= lam if g==1 else (pi**2 if g==2 else zeta(g))
    return v
cw = {int(w): mpf(v) for w, v in json.load(open('cw_dps420.json')).items()}
def lll_relation(vals, D):
    # vals: list of mpf; find integer relation including vals[-1]
    n = len(vals)
    S = mpmath.power(10, D)
    M = []
    for i in range(n):
        row = [0]*n + [int(mpmath.nint(vals[i]*S))]
        row[i] = 1
        M.append(row)
    A = fmpz_mat([[fmpz(x) for x in row] for row in M])
    B = A.lll(delta=mpmath.mpf('0.99'), gram='exact')
    rel = [int(B[0, i]) for i in range(n)]
    resid = sum(mpf(rel[i])*vals[i] for i in range(n))
    return rel, resid, int(B[0, n])
for w in [6]:
    basis = monoms(w)
    vals = [mval(m) for m in basis] + [cw[w]]
    rel, resid, last = lll_relation(vals, 390)
    print('w=',w,'basis=',basis)
    print('rel=', rel)
    print('residual=', mpmath.nstr(resid, 8), 'last-coord=', last)
