import os
base = os.getcwd()

step_02 = '''# step_02.py -- D12c Task 1: exact 1-D representation of L_q(p) + 40-digit numerical checks
# L_q(p) = (4/p) ∫_0^1 (t^{q/p-q} + (1-t)^{q/p-q})^{1/q} dt      [t = |x|^p]
# stable half form: = (8/p) ∫_0^{1/2} t^{1/p-1} (1+(t/(1-t))^{q(1-1/p)})^{1/q} dt
import mpmath as mp, json

def Lq_stable(p, q, dps, deg=10):
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    eps = 1/p; a = q*(1-eps)
    f = lambda t: t**(eps-1)*(1+(t/(1-t))**a)**(1/q)
    return (8/p)*mp.quad(f, [0, mp.mpf('0.25'), mp.mpf('0.5')], maxdegree=deg)

def Lq_xGL(p, q, dps, npanels=100):
    """Independent: direct ell^q arclength, x-parametrisation on half quadrant."""
    mp.mp.dps = dps
    p, q = mp.mpf(p), mp.mpf(q)
    xhalf = mp.mpf('0.5')**(1/p)
    def dsq(x):
        r = (x**p)/(1-x**p)
        return (1 + r**(q*(p-1)/p))**(1/q)
    pts = [xhalf*i/npanels for i in range(npanels+1)]
    return 8*mp.quad(dsq, pts, method='gauss-legendre', maxdegree=10)

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

if __name__ == '__main__':
    mp.mp.dps = 60
    vals = {}
    for (p,q) in [(3,2),(3,3),(3,1.5),(5,4)]:
        a, b = Lq_stable(p,q,220), Lq_xGL(p,q,220)
        vals[f"{p},{q}"] = dict(value=mp.nstr(a,70), diff=mp.nstr(abs(a-b),3))
        print(f"(p,q)=({p},{q}): {mp.nstr(a,60)}  |diff|={mp.nstr(abs(a-b),3)}")
    print("L_2(3)=P(3):", mp.nstr(abs(Lq_stable(3,2,220)-P_squircle(3,220)),3))
    print("L_2(5)=P(5):", mp.nstr(abs(Lq_stable(5,2,220)-P_squircle(5,220)),3))
    print("L_3(3)=2pi_3:", mp.nstr(abs(Lq_stable(3,3,220)-2*pip_cosh(3,220)),3))
    print("L_5(5)=2pi_5:", mp.nstr(abs(Lq_stable(5,5,220)-2*pip_cosh(5,220)),3))
    print("L_1(3)=8:", mp.nstr(abs(Lq_stable(3,1,80)-8),3))
    json.dump(vals, open('task1_values.json','w'), indent=1)
'''
open(os.path.join(base,'step_02.py'),'w').write(step_02)

step_03 = '''# step_03.py -- D12c: u-substitution form, exact split L_q(p)/8 = 2^{-1/p} + eps*H_q(eps),
# and the sup-norm ruler: L_infty(p) = 8*2^{-1/p} exactly.
import mpmath as mp

def Gq_split(eps, q, dps=80):
    """L_q(1/eps)/8 = 2^{-eps} + eps ∫_0^1 u^{eps-1}(1+u)^{-1-eps}[(1+u^{q(1-eps)})^{1/q}-1] du"""
    mp.mp.dps = dps
    eps, q = mp.mpf(eps), mp.mpf(q)
    g = lambda u: u**(eps-1)*(1+u)**(-1-eps)*((1+u**(q*(1-eps)))**(1/q)-1)
    return 2**(-eps) + eps*mp.quad(g, [0, mp.mpf('0.5'), 1], maxdegree=10)

if __name__ == '__main__':
    import json
    mp.mp.dps = 60
    vals = json.load(open('task1_values.json'))
    for key, q in [("3,2",2),("3,3",3),("3,1.5",1.5),("5,4",4)]:
        p = mp.mpf(key.split(',')[0]); known = mp.mpf(vals[key]['value'])
        print(key, "|split - known| =", mp.nstr(abs(Gq_split(1/p,q)*8-known),3))
    # sup-norm ruler: L_inf(p) = 8*2^{-1/p}; two derivations agree
    # (i) geometric: quarter arc has max=|dx| on [0,1/2], gives 8*2^{-1/p}
    # (ii) q->inf of the u-form: 8*eps∫_0^1 u^{eps-1}(1+u)^{-1-eps}du = 8*2F1(1+eps,eps;1+eps;-1) = 8*2^{-eps}
    for p in [3,5,10]:
        e = 1/mp.mpf(p)
        val = 8*e*mp.quad(lambda u: u**(e-1)*(1+u)**(-1-e), [0,1], maxdegree=10)
        print(f"L_inf({p}): u-form {mp.nstr(val,30)} vs 8*2^-1/p {mp.nstr(8*2**(-e),30)} diff {mp.nstr(abs(val-8*2**(-e)),3)}")
'''
open(os.path.join(base,'step_03.py'),'w').write(step_03)
print("wrote step_02.py, step_03.py")
