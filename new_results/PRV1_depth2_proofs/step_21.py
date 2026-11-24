
import mpmath as mp
import sympy as sp
mp.mp.dps = 100

def depth2(x, y, p=2):
    rho = max(abs(x), abs(x*y)); assert 0 < rho < 1
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

ts = [mp.mpf('0.13'), mp.mpf('0.29'), mp.mpf('0.41'), mp.mpf('0.57'), mp.mpf('0.83')]
Z2, Z3 = mp.zeta(2), mp.zeta(3)
M = mp.matrix([[Z2*mp.log(t), Z2*mp.log(1-t), Z2*mp.log(1+t), mp.mpf(1)] for t in ts[:4]])
bvec = mp.matrix([H_dec(t)-E2(t) for t in ts[:4]])
sol = mp.lu_solve(M, bvec)
print("coefficients:")
for name, v in zip(['a (Z2 log t)','b (Z2 log(1-t))','c (Z2 log(1+t))','const'], sol):
    print('  ', name, '=', mp.nstr(v, 30))
# verify at 5th point
t5 = ts[4]
pred = sol[0]*Z2*mp.log(t5)+sol[1]*Z2*mp.log(1-t5)+sol[2]*Z2*mp.log(1+t5)+sol[3]
print("5th point residual:", mp.nstr((H_dec(t5)-E2(t5))-pred, 15))
# nsimplify
for name, v in zip(['a','b','c'], sol[:3]):
    print(name, '=', sp.nsimplify(mp.nstr(v, 60)))
print("const in basis {Z2 log2, Z3}:")
const = sol[3]
# solve c1*Z2*log2 + c2*Z3 = const
det = Z2*mp.log(2)*Z3
c1 = const/(Z2*mp.log(2))  # try pure Z2 log2
print('  const/(Z2 log2) =', mp.nstr(c1, 30))
print('  const/Z3 =', mp.nstr(const/Z3, 30))
print('  nsimplify vs pi^2 log 2 and zeta(3):', sp.nsimplify(mp.nstr(const, 60), [sp.pi**2*sp.log(2), sp.zeta(3)]))
