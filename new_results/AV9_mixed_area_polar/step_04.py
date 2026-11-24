
# step_04.py : NUMERICAL verification of exact c_1..c_8 (two independent methods)
import sympy as sp
import mpmath as mp
import av9core as C
mp.mp.dps = 110

N = 8
h, e = C.selfdual_coeffs(N)
tab = C.beta_int_table(N)
c_exact = {n: C.integrate_poly(4*e[n], tab) for n in range(1, N+1)}

# (a) numerical integration of coefficient functions
lv, lw = C.lv, C.lw
cnum_a = {}
for n in range(1, N+1):
    f = sp.lambdify((lv, lw), e[n], 'mpmath')
    F = lambda v: 4*f(mp.log(v), mp.log(1-v))
    cnum_a[n] = mp.quad(F, [0, mp.mpf('0.5'), 1])
    d = abs(cnum_a[n] - mp.mpf(sp.N(c_exact[n], 100)))
    print(f"n={n}: exact-vs-quad diff = {mp.nstr(d,3)}")

# (b) Cauchy-FFT from V(eta) samples
def Vnum(etaval):
    q = 1+etaval
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [mp.mpc(0), mp.mpc('0.5'), mp.mpc(1)])
M = 24
rho = mp.mpf('0.02')
samples = []
for k in range(M):
    z = rho*mp.exp(2*mp.pi*mp.j*k/M)
    samples.append(Vnum(z))
cnum_b = []
for n in range(N+1):
    s = sum(samples[k]*mp.exp(-2*mp.pi*mp.j*k*n/M) for k in range(M))/M/rho**n
    cnum_b.append(s)
print("c_0 (should be 4):", mp.nstr(cnum_b[0], 20))
for n in range(1, N+1):
    d = abs(cnum_b[n] - mp.mpf(sp.N(c_exact[n], 100)))
    print(f"n={n}: FFT-vs-exact diff = {mp.nstr(d,3)}   value={mp.nstr(mp.re(cnum_b[n]),15)}")
