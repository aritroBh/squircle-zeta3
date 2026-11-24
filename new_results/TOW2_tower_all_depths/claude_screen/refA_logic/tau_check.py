"""Independent exact check of tau1_n, tau2_n (n<=7) as Taylor coefficients of DL1, DL2 at v=i, using sympy series (not the closed formula)."""
import sympy as sp, sys
sys.path.insert(0,'.')
import tow2lib as T
v = sp.symbols('v'); i = sp.I
a = 1+sp.sqrt(2); b = sp.sqrt(2)-1
DL1 = 1/(v-b) + 1/(v+a) - 1/v
DL2 = 1/(v-a) + 1/(v+b) - 1/v
ok = True
for n in range(0, 8):
    # coefficient of (v-i)^{n+1}
    c1 = sp.diff(DL1, v, n+1).subs(v, i)/sp.factorial(n+1)
    c2 = sp.diff(DL2, v, n+1).subs(v, i)/sp.factorial(n+1)
    d1 = sp.nsimplify(sp.simplify(sp.radsimp(sp.expand(c1)) - T.tau1f(n)))
    d2 = sp.nsimplify(sp.simplify(sp.radsimp(sp.expand(c2)) - T.tau2f(n)))
    z1 = sp.simplify(sp.expand(sp.radsimp(c1 - T.tau1f(n)))) == 0
    z2 = sp.simplify(sp.expand(sp.radsimp(c2 - T.tau2f(n)))) == 0
    print(n, 'tau1', z1, 'tau2', z2, '| conj(tau1_n) == (-1)^n tau2_n:', sp.simplify(sp.conjugate(T.tau1f(n)) - (-1)**n*T.tau2f(n)) == 0)
    ok &= z1 and z2
# DL1(i)=DL2(i)=0
print('DL1(i), DL2(i):', sp.simplify(DL1.subs(v,i)), sp.simplify(DL2.subs(v,i)))
print('ALL tau OK' if ok else 'TAU MISMATCH')
