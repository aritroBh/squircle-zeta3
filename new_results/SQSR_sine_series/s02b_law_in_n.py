# Rewrite the exact b_n'' law of s02 in terms of n, with
#  class n=4t+1:  h = H_{(n-1)/4}, g = H_{(n-1)/2}, tau = sum_{i>(n-1)/4} 1/(2i-1)^2
#  class n=4t+3:  h = H_{(n+1)/4}, g = H_{(n+1)/2}, tau = sum_{i>(n+1)/4} 1/(2i-1)^2
import pickle, sympy as sp
from symsum import t, LN2, ZETA, Hsym, H_shift
n, h, g, tau, L2, PI = sp.symbols('n h g tau ln2 pi')
law = pickle.load(open('b2_law.pkl', 'rb'))
for e in (1, 3):
    ex = law[e]
    if e == 1:
        sub = {Hsym(1, 1): h, Hsym(1, 2): g,
               Hsym(2, 2): (Hsym(2, 1)/4 + (PI**2/8 - tau)),   # H2_{2t} = H2_t/4 + O2_t, O2_t = pi^2/8 - tau
               }
        ts = (n-1)/4
    else:
        # H_{t} = H_{t+1} - 1/(t+1), H_{2t} = H_{2t+2} - 1/(2t+1) - 1/(2t+2); O2_{t+1} = pi^2/8 - tau
        sub = {Hsym(1, 1): h - 1/(t+1), Hsym(1, 2): g - 1/(2*t+1) - 1/(2*t+2),
               Hsym(2, 2): (Hsym(2, 1) + 1/(t+1)**2)/4 + (PI**2/8 - tau) - 1/(2*t+2)**2 - 1/(2*t+1)**2 + 0}
        ts = (n-3)/4
    ex2 = sp.expand(ex.subs(sub).subs({ZETA[2]: PI**2/6, LN2: L2}))
    ex2 = sp.expand(ex2.subs(t, ts))
    coll = sp.collect(ex2, [L2, PI, tau, h, g, Hsym(2, 1)], evaluate=False)
    print('n = %d mod 4:' % e)
    for k, v in coll.items():
        print('   [%s]  %s' % (k, sp.factor(sp.simplify(v))))
