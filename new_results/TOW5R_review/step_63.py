code = '''"""step_07b (TOW5R): from-scratch L_r = (log A)^{(r)}(i)/r! at r = 9, 10, integer keys.
NO tow2lib/tow5lib formulas: my own S_m(t) (TOW1 def), my own exact Q(i) series arithmetic.
  log A = (2d-1) log v + log S_m(t(v)) + (2-2d) log(v+i) - log(v^2-a^2) - log(v^2-b^2),
  t(v) = v/4 + 1/2 - 1/(4v)   [exact series about v = i to any order].
Series in u = v - i over Q(i) (Fraction pairs). Compare vs tow5lib._L[r] and tow2lib.L(r) evaluated.
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE + '/tow5_copy')
import tow2lib as T2
import tow5lib as T5
d, m = T2.d, T2.m
t0 = time.time()

# ---- exact Q(i) series arithmetic: element = list of (re, im) Fractions ----
G0 = (Fr(0), Fr(0)); G1 = (Fr(1), Fr(0)); GI = (Fr(0), Fr(1))
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def ginv(a):
    dn = a[0]*a[0]+a[1]*a[1]
    return (a[0]/dn, -a[1]/dn)
def sadd(A, B):
    n = max(len(A), len(B)); A = A + [G0]*(n-len(A)); B = B + [G0]*(n-len(B))
    return [gadd(a, b) for a, b in zip(A, B)]
def smul(A, B, D):
    n = min(len(A)+len(B)-1, D+1)
    out = [G0]*n
    for a_, ca in enumerate(A):
        if a_ > D: break
        if ca == G0: continue
        for b_ in range(min(len(B), n-a_)):
            if B[b_] != G0: out[a_+b_] = gadd(out[a_+b_], gmul(ca, B[b_]))
    return out
def sscale(A, c): return [gmul(a, c) for a in A]
def sdiv(A, B, D):
    out = [G0]*(D+1)
    for n in range(D+1):
        acc = A[n] if n < len(A) else G0
        for k in range(1, n+1):
            if k < len(B) and k <= n and out[n-k] != G0 and B[k] != G0:
                acc = gsub(acc, gmul(B[k], out[n-k]))
        out[n] = gmul(acc, ginv(B[0]))
    return out
def slog(A, D):
    # log(A), A[0] != 0; result has log(A[0]) symbolic -- we only need derivatives, so compute
    # logA' = A'/A series to D-1, integrate (constant term irrelevant)
    Ap = [sscale(A[k+1:k+2], (Fr(k+1),Fr(0)))[0] if k+1 < len(A) else G0 for k in range(D)]
    q = sdiv(Ap, A, D-1)
    out = [G0]*(D+1)
    for k in range(1, D+1):
        out[k] = gmul(q[k-1], (Fr(1, k), Fr(0)))
    return out

# ---- S_m(t) as exact polynomial in t over Q: my own (from step_05) ----
t = sp.Symbol('t')
def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-jj) * comb(k, jj) * jj**n for jj in range(k+1)) // factorial(k)
def wfun(s, mm):
    if s < 0 or s > mm: return Fr(0)
    return Fr(factorial(s) * stir2(mm, s) * 2**mm, factorial(mm))
def S_poly_own(dd, mm):
    q2 = t**2 + (1-t)**2
    out = sp.Integer(0)
    for j in range(1, dd+1):
        cj = sp.binomial(sp.Rational(1,2), j)
        inner = sp.Integer(0)
        for s in range(0, j+1):
            w1 = wfun(s, mm); w2 = wfun(j-s, dd-mm)
            if w1 == 0 or w2 == 0: continue
            inner += sp.binomial(j, s) * sp.Rational(w1) * sp.Rational(w2) * t**(2*(j-s)) * (1-t)**(2*s)
        out += cj * q2**(dd-j) * inner
    return sp.Poly(sp.expand(out), t)

# ---- series of t(v) about v = i: t(v) = v/4 + 1/2 - 1/(4v) ----
DMAX = 11
Q4 = (Fr(1,4), Fr(0))
vser = [ (Fr(0),Fr(1)), G1 ] + [G0]*(DMAX-1)          # v = i + u
vinv = [GI]*(1)  # placeholder
# 1/v = 1/(i+u): series
onev = [G0]*(DMAX+1); onev[0] = ginv(GI)
cur = G1
# 1/(i+u) = (1/i) * 1/(1 + u/i) = sum (-1)^n u^n / i^{n+1}
for n in range(DMAX+1):
    # coefficient = (-1)^n * i^{-(n+1)}
    ip = [GI, G0]  # use gmul loop
# simpler: use sympy-free power: i^n cycles
def gpow_i(n):
    n %= 4
    return [G1, GI, (Fr(-1),Fr(0)), (Fr(0),Fr(-1))][n]
onev = [gmul(gpow_i(-(n+1)), (Fr((-1)**n),Fr(0))) for n in range(DMAX+1)]
tser = [G0]*(DMAX+1)
for n in range(DMAX+1):
    c = Fr(0)
    if n == 1: c += Fr(1,4)
    if n == 0: c += Fr(1,2)
    tser[n] = (c, Fr(0))
tser = [gsub(tser[n], gmul(onev[n], Q4)) for n in range(DMAX+1)]
tstar0 = (Fr(1,2), Fr(1,2))
assert tser[0] == tstar0 and tser[1] == G0, (tser[0], tser[1])   # t(i) = t*, t'(i) = 0

def Lr_own(dd, mm, r):
    D = r + 1
    SP = S_poly_own(dd, mm)
    coeffs = [Fr(sp.Rational(SP.nth(j))) for j in range(SP.degree()+1)]
    # S(t(v)): compose poly with tser
    Sv = [G0]*(D+1)
    tp = [G1] + [G0]*D
    for j, cj in enumerate(coeffs):
        if cj != 0:
            Sv = sadd(Sv, sscale(tp, (cj, Fr(0))))
        tp = smul(tp, tser, D)
    logS = slog(Sv, D)
    # algebraic parts as series of log-derivatives: build logA directly from pieces
    # (2d-1) log v: series of log(v) about i: log(i+u) = log i + sum (-1)^{n+1} u^n/(n i^n)
    logv = [G0]*(D+1)
    for n in range(1, D+1):
        logv[n] = gmul(gpow_i(-n), (Fr((-1)**(n+1), n), Fr(0)))
    # (2-2d) log(v+i): log(2i+u)
    logvi = [G0]*(D+1)
    for n in range(1, D+1):
        logvi[n] = gmul(gpow_i(-n), (Fr((-1)**(n+1), n), Fr(0)))
        logvi[n] = gmul(logvi[n], gpow_i(-0))  # same form: 1/(n (2i)^n) with sign: (2i) power
        logvi[n] = gmul((Fr((-1)**(n+1), n), Fr(0)), [gpow_i(-n), G0][0])
        # (v+i) = 2i + u = 2i(1 + u/(2i)): log = log(2i) + sum (-1)^{n+1} u^n/(n (2i)^n)
        logvi[n] = gmul((Fr((-1)**(n+1), n), Fr(0)), gmul(gpow_i(-n), (Fr(1, 2**n), Fr(0))))
    # -log(v^2 - c^2), c in {a, b}: v^2 - c^2 at v=i: -1 - c^2; series of v^2: (i+u)^2
    v2 = smul(vser, vser, D)
    a2 = (Fr(3), Fr(0)); r2c = None
    # work in Q(i,sqrt2) for this piece? AVOID sqrt2: use the pair (v^2-a^2)(v^2-b^2) = v^4-6v^2+1 = G(v)
    # -log(v^2-a^2)-log(v^2-b^2) = -log G(v), G = v^4 - 6 v^2 + 1  -- all in Q(i)!
    v4 = smul(v2, v2, D)
    Gser = [gsub(gadd(v4[n], sscale(v2[n:n+1], (Fr(-6),Fr(0)))[0] if n < len(v2) else G0), (Fr(-1 if n==0 else 0),Fr(0))) for n in range(D+1)]
    logG = slog(Gser, D)
    logA = [G0]*(D+1)
    for n in range(1, D+1):
        term = gadd(sscale(logv[n:n+1], (Fr(2*dd-1),Fr(0)))[0], logS[n])
        term = gadd(term, sscale(logvi[n:n+1], (Fr(2-2*dd),Fr(0)))[0])
        term = gsub(term, logG[n])
        logA[n] = term
    return logA[r]  # = (log A)^{(r)}(i)/r!

import pickle
bad = []; tested = 0
keys = [(dd, mm) for dd in range(7, 15) for mm in range((dd+1)//2, dd+1)] + [(18,18),(18,12),(22,22),(22,13),(26,26),(26,14),(30,30)]
for (dd, mm) in keys:
    for r in [9, 10]:
        own = Lr_own(dd, mm, r)
        ref = T5._L[r].subs({d: dd, m: mm})
        refc = (Fr(sp.re(ref)), Fr(sp.im(ref)))
        ref2 = T2.L(r).subs({d: dd, m: mm})
        ref2c = (Fr(sp.re(sp.expand(ref2))), Fr(sp.im(sp.expand(ref2))))
        tested += 1
        if own != refc or own != ref2c:
            bad.append((dd, mm, r, own, refc))
    print(f'key ({dd},{mm}) r=9,10 ok ({round(time.time()-t0,1)}s)', flush=True)
print(f'STEP_07B VERDICT: from-scratch L_r (own S + own Q(i) series, -log G for the a/b part): {tested} values, bad = {len(bad)}', flush=True)
for b in bad[:4]: print('  BAD', b)
assert not bad
'''
wd = '<project>/research/08_theory/TOW5R_review'
open(wd+'/step_07b.py','w').write(code)
import subprocess, datetime
PY = '<project>/.venv/bin/python'
t0 = datetime.datetime.now()
r = subprocess.run(['nice','-n','10',PY,'step_07b.py'], cwd=wd, capture_output=True, text=True, timeout=1700)
print(r.stdout[-2000:])
print('STDERR tail:', r.stderr[-500:], 'rc:', r.returncode, 'elapsed:', datetime.datetime.now()-t0)
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_07b: from-scratch L(9),L(10) own series vs T5 cache + T2.L rc={r.returncode} ({(datetime.datetime.now()-t0).seconds}s)\n")
