import json, mpmath
mpmath.mp.dps = 405
from mpmath import mpf, log, pi, zeta
from flint import fmpz_mat, fmpz
lam = log(2)
def monoms(w):
    out = []
    def rec(rem, factors, start):
        if rem == 0: out.append(list(factors)); return
        for g in range(start, rem+1):
            if g in (1,2) or g%2==1: rec(rem-g, factors+[g], g)
    rec(w, [], 1)
    return out
def mval(m):
    v = mpf(1)
    for g in m: v *= lam if g==1 else (pi**2 if g==2 else zeta(g))
    return v
cw = {int(w): mpf(v) for w, v in json.load(open('cw_dps420.json')).items()}
def lll_relation(vals, D):
    n = len(vals); S = mpmath.power(10, D)
    M = []
    for i in range(n):
        row = [0]*n + [int(mpmath.nint(vals[i]*S))]; row[i] = 1; M.append(row)
    A = fmpz_mat([[fmpz(x) for x in row] for row in M])
    B = A.lll(delta=mpmath.mpf('0.99'), gram='exact')
    return [int(B[0, i]) for i in range(n)], sum(mpf(int(B[0,i]))*vals[i] for i in range(n))
import math
results = {}
for w in [5,7,8,9]:
    basis = monoms(w)
    vals = [mval(m) for m in basis] + [cw[w]]
    rel, resid = lll_relation(vals, 390)
    g = 0
    for x in rel: g = math.gcd(g, abs(x))
    reln = [x//g for x in rel]
    results[w] = (basis, reln)
    print('w=',w,'dim=',len(vals),'rel=', reln)
    print('   residual=', mpmath.nstr(resid, 8), ' maxcoeff=', max(abs(x) for x in rel))
json.dump({str(w): {'basis': results[w][0], 'rel': results[w][1]} for w in results}, open('lll_relations_dps400.json','w'))
with open('progress.log','a') as f:
    f.write('step_24: LLL identifications at dps 400 (D=390) for w=5..9 -> lll_relations_dps400.json\n')
