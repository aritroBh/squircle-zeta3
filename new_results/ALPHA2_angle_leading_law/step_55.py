import sys, os, pickle, json
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from mpmath import mp, mpf, nstr, pi, acos, gamma, hyper, diffs, factorial
mp.dps = 55
from assemble import beta_table, N_integrals
import sympy as sp
A = [1,1,1,2,5,16,61,272,1385]
def eps(m): return (-1)**((m+1)//2)
fits = pickle.load(open('Fd_fits.pkl','rb'))
t = sp.symbols('t')
results = {}
names = ['c12','c35','c13','c45','c23','c57','c14','c67']
cvals = {'c12':'1/2','c35':'3/5','c13':'1/3','c45':'4/5','c23':'2/3','c57':'5/7','c14':'1/4','c67':'6/7'}
for name in names:
    c, r, tab = beta_table(name)
    N = N_integrals(r, dmax=6, dps=50)
    ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
    alpha = acos(ca)
    def Eafun(s):
        return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
    ew = [list(diffs(Eafun, 0, n))[-1]/factorial(n) for n in range(8)]
    tv = (1-ca)/(1+ca)
    def Fd(d):
        if d == 1: return mpf(1)
        return mpf(str(sp.N(fits[d].subs(t, tv), 55)))
    row = {}
    for w in range(2, 8):
        d = w-1
        tot = mpf(0)
        for (dd, key, rho), beta in tab.items():
            if dd != d: continue
            rhom = mpf(rho.numerator)/rho.denominator
            tot += mpf(beta.numerator)/beta.denominator * N[(key, rhom)]
        trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
        pred = eps(w-2)*A[w-1]*Fd(d)*trig*ew[w]/factorial(w-2)
        row[w] = (nstr(tot,30), nstr(pred,30), nstr(abs(tot-pred),3))
    results[name] = row
    print(f"c={cvals[name]}: " + "  ".join(f"w{w}:|d|={row[w][2]}" for w in range(2,8)))
pickle.dump(results, open('lawcheck.pkl','wb'))
