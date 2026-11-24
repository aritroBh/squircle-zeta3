import sys, os, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from mpmath import mp, mpf, nstr, pi, acos, gamma, hyper, diffs, factorial
mp.dps = 55
from fractions import Fraction as Fr
from math import factorial as mfac
from assemble import N_integrals
import sympy as sp
def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq
def extract_beta18(name):
    data = pickle.load(open(f'kern18_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, 9):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 19):
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                beta = top[0]/mfac(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab
A = [1,1,1,2,5,16,61,272,1385, 7936]
def eps(m): return (-1)**((m+1)//2)
fits = pickle.load(open('Fd_fits.pkl','rb'))
t = sp.symbols('t')
for name in ['c12','c23']:
    c, r, tab = extract_beta18(name)
    N = N_integrals(r, dmax=7, dps=50)
    ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
    def Eafun(s):
        return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
    ew = [list(diffs(Eafun, 0, n))[-1]/factorial(n) for n in range(9)]
    tv = (1-ca)/(1+ca)
    def Fd(d):
        if d == 1: return mpf(1)
        return mpf(str(sp.N(fits[d].subs(t, tv), 55)))
    for w in [8]:
        d = w-1
        tot = mpf(0)
        for (dd, key, rho), beta in tab.items():
            if dd != d: continue
            rhom = mpf(rho.numerator)/rho.denominator
            tot += mpf(beta.numerator)/beta.denominator * N[(key, rhom)]
        trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
        pred = eps(w-2)*A[w-1]*Fd(d)*trig*ew[w]/factorial(w-2)
        print(f"{name} w={w}: asm={nstr(tot,25)} pred={nstr(pred,25)} |diff|={nstr(abs(tot-pred),3)}")
