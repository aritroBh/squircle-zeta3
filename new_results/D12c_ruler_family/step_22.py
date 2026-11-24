import mpmath as mp

def Lq_stable(p, q, dps, deg=10):
    """Stable: L_q(p) = (8/p) ∫_0^{1/2} t^{eps-1} (1+(t/(1-t))^{q(1-eps)})^{1/q} dt"""
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    eps = 1/p
    a = q*(1-eps)
    f = lambda t: t**(eps-1)*(1+(t/(1-t))**a)**(1/q)
    return (8/p)*mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5')], maxdegree=deg)

def Lq_xGL(p, q, dps, npanels=100):
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    xhalf = mp.mpf('0.5')**(1/p)
    def dsq(x):
        r = (x**p)/(1-x**p)
        return (1 + r**(q*(p-1)/p))**(1/q)
    pts = [xhalf*i/npanels for i in range(npanels+1)]
    return 8*mp.quad(dsq, pts, method='gauss-legendre', maxdegree=10)

import json
mp.mp.dps = 60
vals = {}
for (p,q) in [(3,2),(3,3),(3,1.5),(5,4)]:
    a = Lq_stable(p,q,220)
    b = Lq_xGL(p,q,220)
    d = abs(a-b)
    vals[f"{p},{q}"] = dict(value=mp.nstr(a,70), diff=mp.nstr(d,3), digits=float(-mp.log10(d)))
    print(f"(p,q)=({p},{q}): {mp.nstr(a,60)}")
    print(f"          |diff| = {mp.nstr(d,3)} ({mp.nstr(-mp.log10(d),5)} digits)")

# Also L_2 = P and L_p = 2*pi_p at high precision with the stable form
def P_squircle(p, dps=220):
    mp.mp.dps = dps
    p = mp.mpf(p)
    def ds(th):
        dx = -(2/p)*mp.cos(th)**(2/p-1)*mp.sin(th)
        dy = (2/p)*mp.sin(th)**(2/p-1)*mp.cos(th)
        return mp.sqrt(dx*dx+dy*dy)
    return 4*mp.quad(ds, [0, mp.pi/16, mp.pi/8, mp.pi/4, 3*mp.pi/8, 7*mp.pi/16, mp.pi/2], maxdegree=10)

def pip_cosh(p, dps=220):
    mp.mp.dps = dps
    h = 1/mp.mpf(p)
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, 1, 4, 16, 64, mp.inf], maxdegree=10)

print()
print("L_2(3) vs P(3): |diff| =", mp.nstr(abs(Lq_stable(3,2,220)-P_squircle(3,220)),3))
print("L_2(5) vs P(5): |diff| =", mp.nstr(abs(Lq_stable(5,2,220)-P_squircle(5,220)),3))
print("L_3(3) vs 2pi_3: |diff| =", mp.nstr(abs(Lq_stable(3,3,220)-2*pip_cosh(3,220)),3))
print("L_5(5) vs 2pi_5: |diff| =", mp.nstr(abs(Lq_stable(5,5,220)-2*pip_cosh(5,220)),3))
json.dump(vals, open('task1_values.json','w'), indent=1)
print("saved task1_values.json")
