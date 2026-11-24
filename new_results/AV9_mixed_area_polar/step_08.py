
# step_08.py : numeric verification of circle-end a_2..a_8 via Cauchy-FFT on V(s)
import mpmath as mp
mp.mp.dps = 100
def Vs(sval):
    q = 2/(1-2*sval)
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [mp.mpc(0), mp.mpc('0.5'), mp.mpc(1)])
exact = {
 2: 3*mp.pi*(mp.pi**2-8)/4,
 4: mp.pi*(15*mp.pi**4 - 352*mp.pi**2 + 1664)/192,
 6: mp.pi*(9699*mp.pi**6 - 96280*mp.pi**4 + 13440*mp.pi**2 - 39936)/11520,
 8: mp.pi*(-418761675*mp.pi**8 + 4482752384*mp.pi**6 - 3451302400*mp.pi**4 - 63651840*mp.pi**2 + 138117120)/77414400,
}
M = 32
rho = mp.mpf('0.005')
samples = [Vs(rho*mp.exp(2*mp.pi*mp.j*k/M)) for k in range(M)]
for n in range(0, 12):
    an = sum(samples[k]*mp.exp(-2*mp.pi*mp.j*k*n/M) for k in range(M))/M/rho**n
    if n in exact:
        d = abs(an - exact[n])
        print(f"a_{n}: FFT={mp.nstr(mp.re(an),25)}  exact={mp.nstr(exact[n],25)}  diff={mp.nstr(d,3)}")
    else:
        print(f"a_{n} (FFT): {mp.nstr(an,8)}")
