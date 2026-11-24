import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)

script = '''"""prove_S34.py -- PROOF for all (d,m) of
   S^{(l)}(t*)/S(t*) = N_l(d,m) / prod_{j=1..l}(2d-2j-1),  l = 3, 4.
Same method as prove_S12.py: normalized derivative = finite sum of explicit rational
functions of (d,m) over Q(i); Stirling tails P_r, r <= 4, proved by recurrence (below).
"""
import sympy as sp
d, m = sp.symbols('d m')
i = sp.I

def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r

def tf(A):
    return (i/2)**(A//2) if A % 2 == 0 else (1+i)/2 * (i/2)**((A-1)//2)
def tf1(A):
    return (-i/2)**(A//2) if A % 2 == 0 else (1-i)/2 * (-i/2)**((A-1)//2)

x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
# recurrence proof for all r <= 4 (symbolic identities)
for r in range(1,5):
    assert sp.expand(P[r].subs(x,x+1) - ((x+1-r)*P[r-1] + P[r])) == 0
    for a in range(0, r+1):
        pass
# anchors: P_r vanishes at x = 0..r (S(r,0)=0 for r>=1 and S(m,m-r)=0 for m<r)
for r in range(1,5):
    for a in range(0, r+1):
        assert P[r].subs(x,a) == 0

def normalized_summand(l, n, r, a, c):
    b = l - a
    jj = d - n; s = m - r
    t = sp.Integer(1)
    for rr in range(0, n):
        R = d - n + 1 + rr
        t *= -2*R/(2*R-3)
    t *= ff(m, r) * ff(d-m, n-r) / ff(d, n)
    t *= P[r].subs(x, m) / ff(m, r)
    t *= P[n-r].subs(x, d-m) / ff(d-m, n-r)
    t *= sp.binomial(l, a) * sp.Rational(2)**n * sp.binomial(n, a-n) * i**(2*n-a) * sp.factorial(a)
    A1 = 2*(r-n) - c; A2 = -2*r - b + c
    t *= sp.binomial(b, c) * ff(2*(jj-s), c) * (-1)**(b-c) * ff(2*s, b-c) * tf(A1) * tf1(A2)
    return sp.expand(t)

def total(l):
    tot = sp.Integer(0)
    for n in range(0, l+1):
        for r in range(0, n+1):
            for a in range(n, l+1):
                for c in range(0, l-a+1):
                    tot += normalized_summand(l, n, r, a, c)
    return sp.expand(tot)

N = {3: (-688)*m + (-672*i)*m**2 + 64*m**3 + (344-280*i)*d + (2768+672*i)*d*m + (-96+864*i)*d*m**2 + (-1384+1144*i)*d**2 + (-3264-864*i)*d**2*m + (-192*i)*d**2*m**2 + (1648-1984*i)*d**3 + (1408+192*i)*d**3*m + (-704+1600*i)*d**4 + (-192)*d**4*m + (96-544*i)*d**5 + (64*i)*d**6,
     4: 13824*i*m - 30784*m**2 - 6912*i*m**3 + 256*m**4 + (-sp.Rational(26576,5)-6912*i)*d + (30784-88320*i)*d*m + (72448+10368*i)*d*m**2 + (-512+6144*i)*d*m**3 + (sp.Rational(60592,3)+44160*i)*d**2 + (-72448+169856*i)*d**2*m + (-54528-9216*i)*d**2*m**2 + (-1024*i)*d**2*m**3 + (-48000-86656*i)*d**3 + (54784-148736*i)*d**3*m + (15872+1536*i)*d**3*m**2 + (sp.Rational(202624,3)+75904*i)*d**4 + (-15872+64768*i)*d**4*m + (-1536)*d**4*m**2 + (-sp.Rational(251904,5)-32640*i)*d**5 + (1536-13312*i)*d**5*m + (sp.Rational(57856,3)+6656*i)*d**6 + (1024*i)*d**6*m + (-3584-512*i)*d**7 + 256*d**8}
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out

import sys, time
l = int(sys.argv[1])
t0 = time.time()
T = total(l)
print(f'l={l}: total built ({round(time.time()-t0,1)}s)', flush=True)
rem = sp.together(T - N[l]/Pi(l))
num = sp.expand(sp.fraction(rem)[0])
print(f'S^{l}/S identity: numerator == 0 : {num==0}  ({round(time.time()-t0,1)}s)', flush=True)
assert num == 0
'''
open('prove_S34.py','w').write(script)
import subprocess, sys
r = subprocess.run(['nice','-n','10',sys.executable,'prove_S34.py','3'], capture_output=True, text=True, timeout=2800)
print(r.stdout.strip(), r.stderr[-150:] if r.returncode else '', flush=True)
