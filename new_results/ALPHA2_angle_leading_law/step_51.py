import sys, os, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from mpmath import mp, mpf, nstr, pi, acos, gamma, hyper, diffs, factorial
mp.dps = 60
from assemble import beta_table, N_integrals
import sympy as sp
c, r, tab = beta_table('c12')
N = N_integrals(r, dmax=6, dps=55)
def assemble_nofac(w):
    d = w-1
    tot = mpf(0)
    for (dd, key, rho), beta in tab.items():
        if dd != d: continue
        rhom = mpf(rho.numerator)/rho.denominator
        tot += mpf(beta.numerator)/beta.denominator * N[(key, rhom)]
    return tot
A = [1,1,1,2,5,16,61,272,1385]
def eps(m): return (-1)**((m+1)//2)
ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
def Eafun(s):
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
ew = [list(diffs(Eafun, 0, n))[-1]/factorial(n) for n in range(9)]
fits = pickle.load(open('Fd_fits.pkl','rb'))
t = sp.symbols('t'); tv = (1-ca)/(1+ca)
def Fd(d):
    if d == 1: return mpf(1)
    return mpf(str(sp.N(fits[d].subs(t, tv), 60)))
print("w : c(w;c=1/2)=Sum beta N  vs  eps A F [trig] e_w/(w-2)!")
for w in range(2, 8):
    d = w-1
    trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
    asm = assemble_nofac(w)
    pred = eps(w-2)*A[w-1]*Fd(d)*trig*ew[w]/factorial(w-2)
    print(f"w={w}: asm={nstr(asm,22)} pred={nstr(pred,22)} |diff|={nstr(abs(asm-pred),3)}")
