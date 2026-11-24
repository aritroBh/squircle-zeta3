import sys, os, pickle
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from mpmath import mp, mpf, nstr, pi, acos, gamma, hyper, diffs, factorial
mp.dps = 60
from assemble import assemble_cw
import sympy as sp
cw, c, r = assemble_cw('c12', wmax=7)
A = [1,1,1,2,5,16,61,272,1385]
def eps(m): return (-1)**((m+1)//2)
ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
alpha = acos(ca)
def Eafun(s):
    c2x = (1+ca)/2
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], c2x)/(1-ca)
ew = [list(diffs(Eafun, 0, n))[-1]/factorial(n) for n in range(9)]
fits = pickle.load(open('Fd_fits.pkl','rb'))
t = sp.symbols('t'); tv = (1-ca)/(1+ca)
def Fd(d):
    if d == 1: return mpf(1)
    return mpf(str(sp.N(fits[d].subs(t, tv), 60)))
for w in range(2, 8):
    d = w-1
    trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
    pred = eps(w-2)*A[w-1]*Fd(d)*trig*ew[w]/factorial(w-2)
    print(f"w={w}: asm={nstr(cw[w],22)} pred={nstr(pred,22)} |diff|={nstr(abs(cw[w]-pred),3)}")
