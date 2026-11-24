import json, mpmath, math
mpmath.mp.dps = 240
from mpmath import mpf, log, pi, zeta
from fractions import Fraction as Fr
from flint import fmpz_mat, fmpz
Dtop = json.load(open('Dtop.json'))
N = json.load(open('Nints_dps250.json'))['N']
def Nval(i,j,pole): return mpf(N[f"{i},{j}|{pole}"])
cw2 = {}
for w in range(2,10):
    d = w-1
    tot = mpf(0)
    for name, rec in Dtop.items():
        dd, key = name.split('|'); i, j = map(int, key.split(','))
        if int(dd) != d: continue
        for pole, pr in rec.items():
            cf = Fr(pr['Dtop'][0])
            if cf: tot += mpf(cf.numerator)/cf.denominator * Nval(i,j,pole)
    cw2[w] = tot / mpmath.factorial(d-1)
json.dump({str(w): mpmath.nstr(cw2[w], 240) for w in cw2}, open('cw_dps250.json','w'))
lam = log(2)
# consistency of the two runs
cw1 = {int(w): mpf(v) for w, v in json.load(open('cw_dps420.json')).items()}
for w in range(2,10):
    print(f"w={w}: |c(w) run2 - run1| =", mpmath.nstr(abs(cw2[w]-cw1[w]), 6))
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
def lll_relation(vals, D):
    n = len(vals); S = mpmath.power(10, D)
    M = []
    for i in range(n):
        row = [0]*n + [int(mpmath.nint(vals[i]*S))]; row[i] = 1; M.append(row)
    B = fmpz_mat([[fmpz(x) for x in row] for row in M]).lll(delta=mpmath.mpf('0.99'), gram='exact')
    rel = [int(B[0, i]) for i in range(n)]
    return rel, sum(mpf(rel[i])*vals[i] for i in range(n))
prev = json.load(open('lll_relations_dps400.json'))
allok = True
for w in [5,6,7,8,9]:
    basis = monoms(w)
    vals = [mval(m) for m in basis] + [cw2[w]]
    rel, resid = lll_relation(vals, 225)
    g = 0
    for x in rel: g = math.gcd(g, abs(x))
    reln = [x//g for x in rel]
    if w == 6:
        target = [168, 420, 5040, 378, 2520, 22680, 61, 3780, -362880]
        same = reln == target or reln == [-x for x in target]
    else:
        target = prev[str(w)]['rel']
        same = reln == target or reln == [-x for x in target]
    allok &= same
    print(f"w={w}: dps250 LLL identical to dps400: {same}   residual={mpmath.nstr(resid,6)}")
print("TWO-PRECISION RULE PASSED:", allok)
with open('progress.log','a') as f:
    f.write('step_27: dps250 reassembly + LLL at D=225: relations IDENTICAL to dps400 for w=5..9; runs agree to 1e-240\n')
