import sympy as sp, pickle, os, sys, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
L0 = sp.log(t); L1t = sp.log(1-t)
def ibp_check(n, typ, verbose=False):
    if typ=='W':
        C = t**(n+1)*w1*w4/(n+1)
        R = -kap1*t**(n+1)*w4/u1 - kap4*t**(n+1)*w1/u4
    elif typ=='M1':
        C = t**(n+1)*L0*w1*w4/(n+1)
        R = -t**n*w1*w4 - kap1*t**(n+1)*L0*w4/u1 - kap4*t**(n+1)*L0*w1/u4
    elif typ=='M2':
        C = (t**(n+1)-1)*L1t*w1*w4/(n+1)
        R = -sum(t**j for j in range(n+1))*w1*w4 - kap1*(t**(n+1)-1)*L1t*w4/u1 - kap4*(t**(n+1)-1)*L1t*w1/u4
    elif typ=='M3':
        C = t**(n+1)*w1**2*w4/(n+1)
        R = -2*kap1*t**(n+1)*w1*w4/u1 - kap4*t**(n+1)*w1**2/u4
    elif typ=='M4':
        C = t**(n+1)*w1*w4**2/(n+1)
        R = -kap1*t**(n+1)*w4**2/u1 - 2*kap4*t**(n+1)*w1*w4/u4
    integrand = {'W': t**n*w1*w4, 'M1': t**n*L0*w1*w4, 'M2': t**n*L1t*w1*w4,
                 'M3': t**n*w1**2*w4, 'M4': t**n*w1*w4**2}[typ]
    lhs = (n+1)*integrand
    rhs = (n+1)*sp.diff(C,t) + R
    d = sp.simplify(sp.expand(lhs - rhs))
    if verbose or d != 0:
        print(typ, n, "residual:", d)
    return d == 0
fails = [(typ,n) for n in range(4) for typ in ['W','M1','M2','M3','M4'] if not ibp_check(n,typ,True)]
print("fails:", fails)
