import sys
sys.path.insert(0, '<project>/research/08_theory/ALPHA2_angle_leading_law')
import os
os.chdir('<project>/research/08_theory/ALPHA2_angle_leading_law')
from assemble import assemble_cw
from mpmath import mp, mpf, nstr, pi, sin, cos, gamma, hyper, diffs, factorial, polylog, sqrt
mp.dps = 60
cw, c, r = assemble_cw('c12', wmax=7)
print("c =", c, "r =", r)
for w, v in cw.items():
    print(f"c({w}; c=1/2) = {nstr(v, 30)}")
# prediction: c(w;a) = eps_{w-2} A_{w-1} F_{w-1} [trig] e_w /(w-2)!
# eps_m = (-1)^floor((m+1)/2); A zigzag; F from fits; trig = (1-ca)/2 even w, (1+ca)/2 odd w
A = [1,1,1,2,5,16,61,272,1385]
def eps(m): return (-1)**((m+1)//2)
ca = 2*mpf(c.numerator)/c.denominator**2*mpf(c.denominator) - 1
ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
alpha = __import__('mpmath').acos(ca)
def Eafun(s):
    c2 = (1+ca)/2
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], c2)/(1-ca)
ew = [list(diffs(Eafun, 0, n))[-1]/factorial(n) for n in range(9)]
import pickle, sympy as sp
fits = pickle.load(open('Fd_fits.pkl','rb'))
t = sp.symbols('t')
tv = (1-ca)/(1+ca)
def Fd(d):
    return mpf(str(sp.N(fits[d].subs(t, tv), 60)))
print("angle check: ca =", nstr(ca, 20), " alpha =", nstr(alpha, 20))
for w in range(2, 8):
    d = w-1
    trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
    pred = eps(w-2)*A[w-1]*Fd(d)*trig*ew[w]/factorial(w-2)
    print(f"w={w}: assembled={nstr(cw[w],25)} pred={nstr(pred,25)} |diff|={nstr(abs(cw[w]-pred),3)}")
