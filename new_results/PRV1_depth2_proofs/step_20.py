
import mpmath as mp
mp.mp.dps = 100

def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y))
    assert 0 < rho < 1
    target = mp.power(10, -mp.mp.dps - 5)
    nmax = int(mp.ceil(mp.log(target*(1-rho))/mp.log(rho)))
    h, xy_power, total = mp.mpf(0), x*y, mp.mpf(0)
    for n in range(1, nmax+1):
        total += h/n**p
        h = x*h + x*xy_power/n
        xy_power *= x*y
    return total

def H_dec(t): return depth2(1-t, -t/(1-t))
def E2(t):
    return (mp.polylog(3, 1-t**2) - mp.polylog(3, 1/(1+t)) + mp.polylog(3, t/(1+t))
            - mp.log(1+t)*mp.polylog(2, -t)/2 - 3*mp.log(1+t)*mp.polylog(2, 1-t)/2
            - mp.log(1+t)*mp.polylog(2, 1-t**2)/4)

# residual R(t) = H_dec - E2 should be a*Z2*log t + b*Z2*log(1-t) + c*Z2*log(1+t) + d*Z2*log2 + e*Z3
ts = [mp.mpf('0.13'), mp.mpf('0.29'), mp.mpf('0.41'), mp.mpf('0.57'), mp.mpf('0.71'), mp.mpf('0.83')]
Z2, Z3 = mp.zeta(2), mp.zeta(3)
rows = []
for t in ts:
    R = H_dec(t) - E2(t)
    rows.append((t, R))
M = mp.matrix([[mp.log(t), mp.log(1-t), mp.log(1+t), mp.log(2), Z3/Z2] for t in ts])
# solve a*Z2*(logs) + e*Z3 = R  -> unknowns (a,b,c,d,e) with basis Z2*log.., Z3
bvec = mp.matrix([R for t, R in rows])
# use first 5 equations
Msq = mp.matrix([[M[i,j] for j in range(5)] for i in range(5)])
bsq = mp.matrix([bvec[i] for i in range(5)])
sol = mp.lu_solve(Msq, bsq)
print("fit coefficients (should be rationals):")
for name, v in zip(['Z2 log t','Z2 log(1-t)','Z2 log(1+t)','Z2 log2','Z3'], sol):
    print(' ', name, '=', mp.nstr(v, 30), ' ~ rational:', mp.nstr(mp.identify(v, ['']), 15) if False else '')
# check residual with 6th point
pred = sum(M[5,j]*sol[j] for j in range(5))
print("6th point residual:", mp.nstr(bvec[5]-pred, 10))
print("coeffs as approx rationals:")
for v in sol:
    print('  ', mp.nstr(v, 25), mp.pslq(mp.matrix([1, v]) if False else [1]))
import sympy as sp
for name, v in zip(['a','b','c','d','e'], sol):
    print(name, '=', sp.nsimplify(mp.nstr(v, 60), rational=True))
