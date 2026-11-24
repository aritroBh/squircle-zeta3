import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/tow3_copy'); os.chdir(os.path.dirname(os.path.abspath(__file__))+'/tow3_copy')
import sympy as sp, tow2lib as T
print('check_tau(n<=7):', T.check_tau(7))
ok = True
for n in range(0, 8):
    c1 = sp.simplify(sp.conjugate(T.tau1(n)) - (-1)**n * T.tau2(n))
    c2 = sp.simplify(sp.conjugate(T.tau2(n)) - (-1)**n * T.tau1(n))
    ok &= (c1==0) and (c2==0)
    if c1!=0 or c2!=0: print('FAIL', n, c1, c2)
print('conj(tau1_n)=(-1)^n tau2_n and conj(tau2_n)=(-1)^n tau1_n for n=0..7:', ok)
