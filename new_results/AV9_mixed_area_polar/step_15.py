
# step_15.py : numeric verification of c_9, c_10 (quad + FFT)
import sympy as sp, mpmath as mp, pickle
import av9core as C
mp.mp.dps = 100
out = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_9_10.pkl','rb')).items()}
N = 10
h, e = C.selfdual_coeffs(N)
for n in [9, 10]:
    f = sp.lambdify((C.lv, C.lw), e[n], 'mpmath')
    val = mp.quad(lambda v: 4*f(mp.log(v), mp.log(1-v)), [0, mp.mpf('0.5'), 1])
    ex = mp.mpf(sp.N(out[n], 90))
    print(f"c_{n}: quad-vs-exact diff = {mp.nstr(abs(val-ex),3)}")
# FFT
def Vnum(etaval):
    q = 1+etaval
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [mp.mpc(0), mp.mpc('0.5'), mp.mpc(1)])
M = 32; rho = mp.mpf('0.015')
samples = [Vnum(rho*mp.exp(2*mp.pi*mp.j*k/M)) for k in range(M)]
for n in [9, 10]:
    cf = sum(samples[k]*mp.exp(-2*mp.pi*mp.j*k*n/M) for k in range(M))/M/rho**n
    ex = mp.mpf(sp.N(out[n], 90))
    print(f"c_{n}: FFT={mp.nstr(mp.re(cf),18)}  exact={mp.nstr(ex,18)}  diff={mp.nstr(abs(cf-ex),3)}")
