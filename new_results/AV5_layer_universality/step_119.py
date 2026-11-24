import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr, matrix, lu_solve
mp.dps = 60
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_hex(c1, c2, p):
    q = p/(p-1)
    def F(z):
        return -spow(c1-z,q-1) - spow(c2-z,q-1) + spow(z,q-1)
    lo, hi = -2*max(abs(c1),abs(c2)), 2*max(abs(c1),abs(c2))
    for _ in range(160):
        mid = (lo+hi)/2
        if F(mid) < 0: lo = mid
        else: hi = mid
    z = (lo+hi)/2
    return (abs(c1-z)**q + abs(c2-z)**q + abs(z)**q)**(1/q)
def h_hexinf(c1,c2):
    if c1*c2 >= 0: return max(abs(c1),abs(c2))
    return abs(c1)+abs(c2)
# bump profile Phi(theta) = (h_p - h_inf)/(q-1), two p values to check convergence
for th in [mpf('0.19634954'), mpf('0.39269908'), mpf('0.58904862')]:
    for p in [mpf(200), mpf(600)]:
        q = p/(p-1)
        Phi = (h_hex(cos(th),sin(th),p) - h_hexinf(cos(th),sin(th)))/(q-1)
        print(f"th={th} p={p}: Phi={nstr(Phi,18)}")
    print()
