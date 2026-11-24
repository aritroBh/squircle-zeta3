import json
from mpmath import mp, mpf, sqrt, sin, cos, pi, nstr, quad, log, exp, diff

mp.dps = 60
path = 'research/08_theory/D15b_polygon_law/vals_m2_8_200.json'
data = json.load(open(path))
for k in sorted(data.keys(), key=int):
    print("key", k, "len", len(data[k]))
    for i, e in enumerate(data[k]):
        print("   ", i, str(e)[:70])

# Direct numerical computation of P(p) and p-derivatives at p=1 for the setup
# K_p = { sum_i |<n_i,x>|^p <= 1 }, n_i = (cos(i pi/4), sin(i pi/4)), i=1..4.
# r(theta) = u(theta)^(-1/p), u = sum |cos(theta - i pi/4)|^p ; P(p) = 8 int_0^{pi/4} sqrt(r^2+r'^2)
# On [0,pi/4]: |c1|=cos(t-pi/4), |c2|=sin(t), |c3|=cos(t+pi/4), |c4|=cos(t)

def make_u(p):
    def u(t):
        return cos(t-pi/4)**p + sin(t)**p + cos(t+pi/4)**p + cos(t)**p
    return u

def P_of_p(p):
    u = make_u(p)
    def integrand(t):
        r = u(t)**(-1/p)
        rp = diff(lambda s: u(s)**(-1/p), t)
        return sqrt(r*r + rp*rp)
    return 8*quad(integrand, [0, pi/4])

P1 = P_of_p(1)
print("P(1) direct =", nstr(P1, 50))
print("octagon per =", nstr(8*(sqrt(2)-1)*sqrt(2-sqrt(2)), 50))

# numerical p-derivatives at p=1 via high-order finite differences (for identification only)
for k in range(1,4):
    dk = diff(P_of_p, 1, n=k, h=mpf('0.01'), method='central')
    print(f"P^{{({k})}}(1) ~", nstr(dk, 25))
