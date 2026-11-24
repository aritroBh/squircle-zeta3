import mpmath as mp
mp.mp.dps = 60

def S_p(theta, p, m):
    A = mp.pi/m
    return mp.fsum([abs(mp.sin(theta + i*A))**p for i in range(m)])
def Sp_p(theta, p, m):
    A = mp.pi/m
    return mp.fsum([p*mp.sin(theta + i*A)**(p-1)*mp.cos(theta + i*A) for i in range(m)])
def P_direct(p, m):
    A = mp.pi/m
    def dens(th):
        S = S_p(th, p, m); Sp = Sp_p(th, p, m)
        r = S**(-mp.mpf(1)/p)
        rp = -S**(-mp.mpf(1)/p - 1)*Sp/p
        return mp.sqrt(r*r+rp*rp)
    return 4*m*mp.quad(dens, [0, A/2], method='tanh-sinh')

# FD of P(p) at p=1 for m=3 with several h, Richardson on h^2 (central) 
m = 3
hs = [mp.mpf('0.01'), mp.mpf('0.005'), mp.mpf('0.0025'), mp.mpf('0.00125')]
D1s = []; D3s = []
for h in hs:
    fp = P_direct(1+h, m); fm = P_direct(1-h, m)
    D1s.append((fp-fm)/(2*h))
    D3s.append((fp - 2*P_direct(1,m) + fm)/h**3)  # not quite P''' ; use proper stencil below
print("central D1 estimates (m=3):", [mp.nstr(v,25) for v in D1s])
# Richardson extrapolate D1 (error ~ h^2, h^4, ...)
def richardson(vals, hs, order0=2):
    # fit polynomial in h^2
    n = len(vals)
    A = mp.matrix([[hs[i]**(2*j) for j in range(n)] for i in range(n)])
    b = mp.matrix(vals)
    c = mp.lu_solve(A, b)
    return c[0]
print("Richardson D1 (m=3):", mp.nstr(richardson(D1s, hs), 45))
print("coef-integral D1: 2.40113226770588735429020690407746529554651722044")

# P''' via 5-point stencil (error h^2): f'''(x) ≈ (-f(x-2h)+2f(x-h)-2f(x+h)+f(x+2h))/(2h^3)
D3s = []
for h in hs:
    f2p = P_direct(1+2*h, m); fp = P_direct(1+h, m)
    fm = P_direct(1-h, m); f2m = P_direct(1-2*h, m)
    D3s.append((-f2m + 2*fm - 2*fp + f2p)/(2*h**3))
print("central D3 estimates:", [mp.nstr(v,20) for v in D3s])
print("Richardson D3:", mp.nstr(richardson(D3s, hs), 40))
print("coef-integral D3: -5.14545501713561583717168113112397507111555356642")
