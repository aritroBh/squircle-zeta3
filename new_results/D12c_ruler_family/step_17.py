import mpmath as mp
mp.mp.dps = 100

def Lq_integral(p, q, dps=100):
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    eps = 1/p
    f = lambda t: (t**(q*(eps-1)) + (1-t)**(q*(eps-1)))**(1/q)
    return 4*eps*mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1], maxdegree=10)

def Lq_arc_x(p, q, dps=100):
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    xhalf = mp.mpf('0.5')**(1/p)
    def dsq(x):
        r = (x**p)/(1-x**p)
        return (1 + r**(q*(p-1)/p))**(1/q)
    quarter = 2*mp.quad(dsq, [0, xhalf/2, xhalf], maxdegree=10)
    return 4*quarter

res = {}
for (p,q) in [(3,2),(3,3),(3,1.5),(5,4)]:
    a = Lq_integral(p,q); b = Lq_arc_x(p,q)
    res[(p,q)] = (a,b,abs(a-b))
    print(f"(p,q)=({p},{q}): |diff| = {mp.nstr(abs(a-b),3)}  (digits agree: {mp.nstr(-mp.log10(abs(a-b)),3)})")

# L_2 = Euclidean perimeter: check at p=3 against independent squircle perimeter formula
# P(p) = 2^(2-1/p) * sqrt(pi) * Gamma(1+1/p) / Gamma(1/2+1/p) * 2F1(...)?  Use direct sqrt-arclength:
def P_squircle(p, dps=100):
    # Euclidean length via theta parametrization with corner-split quadrature
    mp.mp.dps = dps
    p = mp.mpf(p)
    def ds(th):
        dx = -(2/p)*mp.cos(th)**(2/p-1)*mp.sin(th)
        dy = (2/p)*mp.sin(th)**(2/p-1)*mp.cos(th)
        return mp.sqrt(dx*dx+dy*dy)
    return 4*mp.quad(ds, [0, mp.pi/16, mp.pi/8, mp.pi/4, 3*mp.pi/8, 7*mp.pi/16, mp.pi/2], maxdegree=10)

for p in [3,5]:
    a = Lq_integral(p,2); b = P_squircle(p)
    print(f"L_2({p}) vs P({p}): |diff| = {mp.nstr(abs(a-b),3)}")
