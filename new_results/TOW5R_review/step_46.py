code = '''"""step_04 (TOW5R): independent re-derivation of the prove_S5 summand at l = 5.
My own hand derivation of S^{(l)}(t*)/S(t*) from the TOW1 definition
  S_m(t) = sum_j binom(1/2,j) q2^{d-j} sum_s binom(j,s) w(s,m) w(j-s,d-m) t^{2(j-s)}(1-t)^{2s},
  w(s,m) = s! Stir2(m,s) 2^m/m!,  q2 = t^2+(1-t)^2 = 2(t-t*)(t-tbar*), t*=(1+i)/2.
Checks (all symbolic, exact):
  [1] q2 == 2(t-t*)(t-tbar*); t*^2 = i/2; tf/tf1(A) == t*^A, (1-t*)^A for A in -12..12.
  [2] (q2^n)^{(a)}(t*) == 2^n binom(n,a-n) i^{2n-a} a!   (n,a <= 5, direct sympy diff).
  [3] binom(1/2,R-1)/binom(1/2,R) == -2R/(2R-3) (symbolic R).
  [4] P_5(n) == Stirling2(n, n-5) for n = 0..14 (sanity of the proved tail).
  [5] w(m-r,m)/2^m == P_r(m)/m^r_ at integer points (m = 0..14, r <= 5).
  [6] MY total(5) (my own summand implementation) == prove_S5's total(5) symbolically;
      and cancel(mine - N5/(3 Pi_5)) == 0.
"""
import sympy as sp, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy')
d, m = sp.symbols('d m')
i = sp.I
t = sp.Symbol('t')
tstar = (1+i)/2; tbar = (1-i)/2
t0 = time.time()

# [1]
q2 = t**2 + (1-t)**2
assert sp.expand(q2 - 2*(t-tstar)*(t-tbar)) == 0
assert sp.simplify(tstar**2 - i/2) == 0 and sp.simplify(tbar**2 + i/2) == 0
def tf(A):
    return (i/2)**(A//2) if A % 2 == 0 else (1+i)/2 * (i/2)**((A-1)//2)
def tf1(A):
    return (-i/2)**(A//2) if A % 2 == 0 else (1-i)/2 * (-i/2)**((A-1)//2)
for A in range(-12, 13):
    assert sp.simplify(tf(A) - tstar**A) == 0, A
    assert sp.simplify(tf1(A) - tbar**A) == 0, A
print('[1] q2 factorization + tf/tf1 values OK', flush=True)

# [2]
for n in range(0, 6):
    for a in range(0, 6):
        lhs = sp.diff(q2**n, t, a).subs(t, tstar)
        rhs = 2**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
        assert sp.simplify(lhs - rhs) == 0, (n, a)
print('[2] (q2^n)^(a)(t*) formula OK (n,a <= 5, direct diff)', flush=True)

# [3]
R = sp.Symbol('R', integer=True, positive=True)
bh = lambda n: sp.binomial(sp.Rational(1,2), n)
ratio = sp.simplify(bh(R-1)/bh(R))
assert sp.simplify(ratio - (-2*R/(2*R-3))) == 0
print('[3] binom(1/2,R-1)/binom(1,2,R) == -2R/(2R-3) OK', flush=True)

# [4]
x = sp.Symbol('x')
P5 = x*(x-1)*(x-2)*(x-3)*(x-4)**2*(x-5)**2*(3*x**2-23*x+38)/11520
for n in range(0, 15):
    assert sp.expand(P5.subs(x, n) - sp.functions.combinatorial.numbers.stirling(n, n-5)) == 0, n
print('[4] P_5(n) == Stirling2(n, n-5), n = 0..14 OK', flush=True)

# [5]
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760, 5: P5}
from math import comb, factorial
def stir2(n, k):
    if n == k == 0: return 1
    if n <= 0 or k <= 0 or k > n: return 0
    return sum((-1)**(k-j) * comb(k, j) * j**n for j in range(k+1)) // factorial(k)
def ff(x0, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x0-u)
    return r
for mm in range(0, 15):
    for r in range(0, 6):
        w = factorial(mm-r)*stir2(mm, mm-r)*2**mm/factorial(mm) if mm-r >= 0 else 0
        assert sp.Rational(w) - P[r].subs(x, mm)/ff(mm, r) == 0, (mm, r)
print('[5] w(m-r,m)/2^m == P_r(m)/m^r_, m = 0..14 OK', flush=True)

# [6] my own summand, written from my derivation
def my_summand(l, n, r, a, c):
    b = l - a
    jj = d - n; s = m - r
    # (q2^n)^(a)(t*) * binom(l,a)
    q2der = 2**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    # binom(1/2,d-n)/binom(1/2,d)
    cratio = sp.Integer(1)
    for R_ in range(d-n+1, d+1):
        pass  # symbolic range: use product form
    cratio = sp.prod([-2*R_/ (2*R_-3) for R_ in range(d-n+1, d+1)]) if False else None
    cratio = sp.Integer(1)
    for k_ in range(n):
        R_ = d - n + 1 + k_
        cratio *= -2*R_/(2*R_-3)
    # binom(d-n, m-r)/binom(d,m)
    bratio = ff(m, r) * ff(d-m, n-r) / ff(d, n)
    # w w / 2^d
    wratio = P[r].subs(x, m)/ff(m, r) * P[n-r].subs(x, d-m)/ff(d-m, n-r)
    # [t^{2(j-s)} (1-t)^{2s}]^{(b)}(t*) * t*^{-2(d-m)} (1-t*)^{-2m}
    A1 = 2*(jj-s) - c - 2*(d-m); A2 = 2*s - (b-c) - 2*m
    bder = sp.binomial(b, c) * ff(2*(jj-s), c) * (-1)**(b-c) * ff(2*s, b-c) * tf(A1) * tf1(A2)
    return sp.expand(sp.binomial(l, a) * q2der * cratio * bratio * wratio * bder)

def my_total(l):
    tot = sp.Integer(0)
    for n in range(0, l+1):
        for r in range(0, n+1):
            for a in range(n, l+1):
                for c in range(0, l-a+1):
                    tot += my_summand(l, n, r, a, c)
    return sp.expand(tot)

Tm = my_total(5)
print('[6] my total(5) built', round(time.time()-t0,1), 's', flush=True)
# prove_S5's total: rebuild via its own module code path (import by exec of the file's functions)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy')
import importlib.util
spec = importlib.util.spec_from_file_location('prove_S5_copy', os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy/prove_S5.py')
# do NOT exec (it would rerun everything); instead re-implement their summand directly:
def their_summand(l, n, r, a, c):
    b = l - a
    jj = d - n; s = m - r
    tt = sp.Integer(1)
    for rr in range(0, n):
        R_ = d - n + 1 + rr
        tt *= -2*R_/(2*R_-3)
    tt *= ff(m, r) * ff(d-m, n-r) / ff(d, n)
    tt *= P[r].subs(x, m) / ff(m, r)
    tt *= P[n-r].subs(x, d-m) / ff(d-m, n-r)
    tt *= sp.binomial(l, a) * sp.Rational(2)**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    A1 = 2*(r-n) - c; A2 = -2*r - b + c
    tt *= sp.binomial(b, c) * ff(2*(jj-s), c) * (-1)**(b-c) * ff(2*s, b-c) * tf(A1) * tf1(A2)
    return sp.expand(tt)
diff_terms = 0
for n in range(0, 6):
    for r in range(0, n+1):
        for a in range(n, 6):
            for c in range(0, 6-a):
                if sp.expand(my_summand(5,n,r,a,c) - their_summand(5,n,r,a,c)) != 0:
                    diff_terms += 1
print('[6] summand-wise differences (mine vs theirs), 126 terms:', diff_terms, flush=True)
assert diff_terms == 0
import pickle
N5 = sp.sympify(pickle.load(open(os.path.dirname(os.path.abspath(__file__))+'/tow5_copy/N5_proved.pkl','rb'))['N5'], locals={'d':d,'m':m,'I':sp.I})
Pi5 = sp.Integer(1)
for j in range(1, 6): Pi5 *= (2*d-2*j-1)
res = sp.cancel(Tm - N5/(3*Pi5))
print('[6] cancel(my total(5) - N5/(3 Pi_5)) == 0 :', res == 0, f'({round(time.time()-t0,1)}s)', flush=True)
assert res == 0
print('STEP_04 VERDICT: my independent Leibniz re-derivation == prove_S5 total(5) == N5/(3 Pi_5), exact', flush=True)
'''
wd = '<project>/research/08_theory/TOW5R_review'
open(wd+'/step_04.py','w').write(code)
print('written')
