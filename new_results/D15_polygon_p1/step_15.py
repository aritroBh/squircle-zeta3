import mpmath as mp
mp.mp.dps = 50
m = 2; A = mp.pi/m; x = A/2
S0 = lambda th: mp.sin(th) + mp.cos(th)
D0 = lambda th: mp.sqrt(1 + ((mp.cos(th)-mp.sin(th))/S0(th))**2)/S0(th)
I1 = mp.quad(D0, [0, x], method='tanh-sinh')
I2 = mp.quad(D0, [mp.mpf('1e-30'), x], method='tanh-sinh')
print("I tanh-sinh [0,x]:", mp.nstr(I1, 30))
print("I tanh-sinh [1e-30,x]:", mp.nstr(I2, 30))
print("exact sqrt2/2:", mp.nstr(mp.sqrt(2)/2, 30))
I3 = mp.quad(D0, [0, x], method='gauss-legendre')
print("I gauss-legendre:", mp.nstr(I3, 30))
# and the version where logs are computed but discarded:
def D0_withlogs(th):
    ls0 = mp.log(mp.sin(th)); ls1 = mp.log(mp.cos(th))
    T1 = mp.sin(th)*ls0 + mp.cos(th)*ls1
    return mp.sqrt(1 + ((mp.cos(th)-mp.sin(th))/S0(th))**2)/S0(th) + 0*T1
I4 = mp.quad(D0_withlogs, [0, x], method='tanh-sinh')
print("I with dummy logs:", mp.nstr(I4, 30))
