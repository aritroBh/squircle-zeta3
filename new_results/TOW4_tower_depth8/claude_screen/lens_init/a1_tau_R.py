# A1: (a) tau_n n<=8 : own DL1/DL2 series (core.TAU1/2, exact) vs tow4lib tables, vs mpmath numeric Taylor of the ORIGINAL partial-fraction DL1/DL2 (50 digits)
#     (b) conjugation rule conj(tau1_n) = (-1)^n tau2_n
#     (c) R_l = S^(l)(t*)/S(t*) (l=1..4) from the DEFINITION of S_m (own engine) vs N_l/Pi_l of tow2lib, ALL (d,m), 1<=d<=40, 0<=m<=d
import sys, time
sys.dont_write_bytecode=True
sys.path.insert(0,'<project>/research/08_theory/TOW4_tower_depth8')
from core import *
import sympy as sp, mpmath as mp
import tow2lib as T2, tow4lib as T4
mp.mp.dps = 60
a = 1+mp.sqrt(2); b = mp.sqrt(2)-1
DL1 = lambda v: 1/(v-b)+1/(v+a)-1/v
DL2 = lambda v: 1/(v-a)+1/(v+b)-1/v
c1 = mp.taylor(DL1, mp.mpc(0,1), 10); c2 = mp.taylor(DL2, mp.mpc(0,1), 10)
ok_t = True; ok_n=True; ok_c=True
for n in range(9):
    t1 = T4.tau1f(n); t2 = T4.tau2f(n)
    ex1 = (sp.Rational(t1.as_real_imag()[0]), sp.Rational(t1.as_real_imag()[1]))
    mine1 = (TAU1[n].r, TAU1[n].i); mine2 = (TAU2[n].r, TAU2[n].i)
    eq1 = (fmpq(int(ex1[0].p), int(ex1[0].q)) == mine1[0] and fmpq(int(ex1[1].p), int(ex1[1].q)) == mine1[1])
    r2,i2 = t2.as_real_imag()
    eq2 = (fmpq(int(sp.Rational(r2).p), int(sp.Rational(r2).q)) == mine2[0] and fmpq(int(sp.Rational(i2).p), int(sp.Rational(i2).q)) == mine2[1])
    # numeric: coefficient of (v-i)^(n+1)
    num1 = c1[n+1]; num2 = c2[n+1]
    m1 = mp.mpc(float(0)) + mp.mpf(int(mine1[0].p))/int(mine1[0].q) + 1j*0
    m1 = mp.mpc(mp.mpf(int(mine1[0].p))/int(mine1[0].q), mp.mpf(int(mine1[1].p))/int(mine1[1].q))
    m2 = mp.mpc(mp.mpf(int(mine2[0].p))/int(mine2[0].q), mp.mpf(int(mine2[1].p))/int(mine2[1].q))
    d1 = abs(num1-m1); d2 = abs(num2-m2)
    ok_n &= (d1 < mp.mpf(10)**-45 and d2 < mp.mpf(10)**-45)
    ok_t &= eq1 and eq2
    ok_c &= (TAU1[n].conj() == TAU2[n]*((-1)**n))
    print(f'n={n}: exact==tow4lib tables: {eq1 and eq2}; |numeric Taylor(partial fractions) - mine| = {mp.nstr(d1,3)}, {mp.nstr(d2,3)}; conj rule {TAU1[n].conj()==TAU2[n]*((-1)**n)}')
print('TAU: exact match with tow4lib', ok_t, '| numeric 45-digit match with original DL1/DL2 Taylor', ok_n, '| conj rule', ok_c)
# tau also vs tow2lib DEFINITION formula (sqrt2 expressions) at n=7,8, simplified
for n in (7,8):
    for nm,f,mine in (('tau1',T2.tau1,TAU1[n]),('tau2',T2.tau2,TAU2[n])):
        diff = sp.simplify(sp.expand(f(n) - (sp.Rational(int(mine.r.p),int(mine.r.q)) + sp.I*sp.Rational(int(mine.i.p),int(mine.i.q)))))
        print(f'  {nm}({n}) definition (T2) - mine = {diff}')

# (c) R_l from definition vs N_l/Pi_l
d_, m_ = T2.d, T2.m
Rf = {l: sp.lambdify((d_,m_), T2.R[l], 'sympy') for l in range(1,5)}
bad=0; cells=0; t0=time.time()
for d in range(1, 41):
    for m in range(0, d+1):
        c = S_taylor(d, m, 4)
        for l in range(1,5):
            Rl = c[l]*factorial(l)/c[0]          # S^(l)(t*)/S(t*)
            ex = sp.nsimplify(Rf[l](sp.Integer(d), sp.Integer(m)))
            exg = sp.expand(ex)
            re_, im_ = exg.as_real_imag()
            re_=sp.Rational(re_); im_=sp.Rational(im_)
            cells+=1
            if not (Rl.r == fmpq(int(re_.p),int(re_.q)) and Rl.i == fmpq(int(im_.p),int(im_.q))): bad+=1
print(f'R_l (l=1..4) from the S_m definition vs N_l/Pi_l: all 1<=d<=40, 0<=m<=d: cells={cells} bad={bad} ({time.time()-t0:.0f}s)')
