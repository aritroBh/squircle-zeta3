
import os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
os.chdir(WORK)
with open('progress.log','a') as f:
    f.write("D13s turn1: copied step_38.py, step_10_J192.py, step_18.py, planner_sweep_ratio.py + review docs.\n"
            "B1 line confirmed: step_10_J192.py:168 `ecw = (rt2/pi).sqrt()*(-c_w/2).exp()` (sqrt(sqrt(2)/pi), needs sqrt(2/pi)).\n"
            "B2 line confirmed: step_38.py:165 `A1m = iv(0, (z2/(1+z2)).upper())` (encloses 1-arctan(zz)/zz, NOT (1-arctan(zz)/zz)/y).\n"
            "Review result.md + progress.log read fully; steps 22,23,28-32,35 read for H1/H2/H3.\n")

# Verify B2-fix math: f(x) = 1 - atan(x)/x is nondecreasing on [0,inf)
import sympy as sp
x = sp.symbols('x', positive=True)
f = 1 - sp.atan(x)/x
fp = sp.simplify(sp.diff(f, x))
print("f'(x) =", fp)
# f'(x) = (atan(x) - x/(1+x^2))/x^2 ; numerator g(x), g(0)=0
g = sp.atan(x) - x/(1+x**2)
gp = sp.simplify(sp.diff(g, x))
print("g'(x) =", gp, "  -- manifestly >= 0 on x>=0:", sp.simplify(gp) == 2*x**2/(1+x**2)**2)
print("g(0) =", g.subs(x,0), "=> g >= 0 => f' >= 0.  f(0+) =", sp.limit(f, x, 0))
print("f increasing on [0,inf): PROVED")
# quick numeric sanity
import mpmath as mp
mp.mp.dps = 40
xs = [mp.mpf(10)**k/3 for k in range(-6, 7)]
vals = [1 - mp.atan(t)/t for t in xs]
print("sample f increasing:", all(vals[i] < vals[i+1] for i in range(len(vals)-1)))
