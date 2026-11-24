import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
from mpmath import mp, mpf, log, log1p, exp, cos, sin, pi
mp.dps = 50
exec(open('av2core.py').read())
p = mpf('0.9'); ep = 1/p; a = 2*ep
def f(th, r):
    cs = cos(th)*sin(th)
    c = cs**a
    lQ = ep*log1p(-r**p)
    Q = exp(lQ)
    one_Q = -mp.expm1(lQ)
    val = one_Q + (r*r*c)*Q
    return val, cs, c, lQ, Q, one_Q
# test the first/last ts nodes
nodes = _ts_nodes(0, pi/2, '0.055', 180)
nodesr = _ts_nodes(0, 1, '0.055', 180)
th0, _ = nodes[0]; th1, _ = nodes[-1]
r0, _ = nodesr[0]; r1, _ = nodesr[-1]
print("th0:", mp.nstr(th0, 10), "th1:", mp.nstr(th1, 10))
print("r0:", mp.nstr(r0, 10), "r1:", mp.nstr(r1, 10))
for (th, r) in [(th0, r0), (th0, r1), (th1, r0), (th1, r1)]:
    val, cs, c, lQ, Q, one_Q = f(th, r)
    print("ok:", mp.nstr(val, 10), type(val))
