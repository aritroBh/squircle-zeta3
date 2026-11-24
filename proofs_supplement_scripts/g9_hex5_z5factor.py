# Exact check of the chain of equalities in the proof of Lemma z5 (Theorem 26), including the factor 120.
#  (1) bracket coefficient of z5 in  sum c_{n,a,b,0} Ihat(n,0;a,b):   c5 = 24 * int sum_{a+b=5,a,b>=1} g_ab dt   (Beta rule (b): 24 B_{n,0})
#      bracket coefficient of z2 z3:                                  c23 = -12 * int (g32+g23) dt            (Beta rule (c))
#      with g_ab the coefficients of Psi built exactly as in g9_edgelaw5.py (its own code is executed up to the comparison).
#      Check: 120*c5 = Z_e^(5) and 120*c23 = Z_e^(2,3)  (so the bracket coefficients are Z/120, not Z).
#  (2) the same two numbers read off the table g9_hex5_table.txt (lines x^0 | 1 | z5^1 and z2^1 z3^1): 120 * coefficient - law = 0.
import os, sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'g9_edgelaw5.py')).read().split('d5 = sp.expand(Z5 - law5)')[0]
ns = {}
exec(compile(src, 'g9_edgelaw5_head', 'exec'), ns)
integ, g, law5, law23 = ns['integ'], ns['g'], ns['law5'], ns['law23']
c5 = 24*integ(sum(g.get((a, 5-a), 0) for a in range(1, 5)))     # = (Z5/L)/120 * L ... g_ab carries the factor L; here in units of L
c23 = -12*integ(g.get((3, 2), 0) + g.get((2, 3), 0))
print('(1) 120*c5/L - law5 =', sp.expand(120*c5 - law5), ';  120*c23/L - law23 =', sp.expand(120*c23 - law23))
assert sp.expand(120*c5 - law5) == 0 and sp.expand(120*c23 - law23) == 0
# (2) table
A0, A1, B, kap, q0, V02, L2 = sp.symbols('A0 A1 B kap q0 V02 L2', positive=True)
loc = dict(A0=A0, A1=A1, B=B, kap=kap, q0=q0, V02=V02, L2=L2, sqrt=sp.sqrt)
found = {}
for line in open(os.path.join(here, 'g9_hex5_table.txt')):
    f = [x.strip() for x in line.split(' | ')]
    if len(f) == 4 and f[0] == 'x^0' and f[1] == '1' and f[2] in ('z2^0 z3^0 z5^1', 'z2^1 z3^1 z5^0'):
        found[f[2]] = sp.sympify(f[3], locals=loc)
rel = {A1: 1 - B, A0: 1 - B - B*kap}
L = sp.sqrt(L2); L0 = -q0/L; L1 = L + q0/L; h2 = V02 - q0**2/L2
xi0, xi1 = A0, A1
xb = (xi0*L0 + xi1*L1)/L; xp_, xt = xi0 + xi1, xi0*xi1; LL = L0*L1/L2; m = h2/L2
Z5 = 360*xt*(h2/L)*(20*xb**3 - 32*xb**2 + 15*xb - 2 - 4*(5*xb - 4)*xt*LL + m*(5*xt*xp_ + 8*xp_**2 - 28*xt - 15*(xp_**2 - 3*xt)*xb))
Z23 = 120*xt*(h2/L)*(12*xb**2 - 12*xb + 2 - 12*(5*xb - 4)*xt*LL + m*(15*xt*xp_ - 3*xp_**2 + 24*xt - 45*xt*xb))
d5 = sp.simplify((120*found['z2^0 z3^0 z5^1'] - Z5).subs(rel)); d23 = sp.simplify((120*found['z2^1 z3^1 z5^0'] - Z23).subs(rel))
print('(2) 120 x table[z5] - Z5 =', d5, ';  120 x table[z2 z3] - Z23 =', d23)
assert d5 == 0 and d23 == 0
print('pass')
