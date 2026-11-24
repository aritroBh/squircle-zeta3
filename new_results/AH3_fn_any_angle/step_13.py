
import mpmath as mp
mp.mp.dps = 60

def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u

def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)

def logR(e, th):
    z = mp.sin(th/2)**2
    if abs(e) < mp.mpf('1e-40'):
        return mp.log(2*z)  # limit
    return mp.log((mp.pi*e/mp.sin(mp.pi*e))*Lambda(e, th))

def a_coeffs(th, J=4):
    z = mp.sin(th/2)**2
    f = lambda e: logR(e, th) - mp.log(2*z)
    return [mp.diff(f, mp.mpf(0), j)/mp.factorial(j) for j in range(1, J+1)]

def cumulants_numeric(th, J=4):
    f = lambda t: mp.log(mp.quad(lambda u: S_of_U(u, th)**(-t), [0,1]))
    return [mp.diff(f, mp.mpf(0), j) for j in range(1, J+1)]

print("Cumulants k_j (of log(1/S_theta)) vs (-1)^j j! a_j from paper's a1..a4:")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.pi/6]:
    aa = a_coeffs(th)
    kk = cumulants_numeric(th)
    for j in range(4):
        target = (-1)**(j+1)*mp.factorial(j+1)*aa[j]
        print(f" th={mp.nstr(th,6)} k{j+1}: numeric={mp.nstr(kk[j],16)}  from_a={mp.nstr(target,16)}  diff={mp.nstr(abs(kk[j]-target),2)}")
    print()

# Also verify closed forms: k1 = 1 + w*ell ; k2 = 1 + pi^2/3 - 2 Li2(z) - (w ell)^2 ; k3 = -12 w S12(z) + 2(1+(w ell)^3)
def S12(z):
    ell = mp.log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-z) + ell*mp.polylog(2, 1-z) + mp.mpf('0.5')*ell**2*mp.log(z)

print("Closed-form cumulants vs numeric:")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.pi/6, mp.mpf('0.9')]:
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1c = 1 + u
    k2c = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3c = -12*w*S12(z) + 2*(1 + u**3)
    kk = cumulants_numeric(th, 3)
    print(f" th={mp.nstr(th,6)}: k1 d={mp.nstr(abs(kk[0]-k1c),2)}  k2 d={mp.nstr(abs(kk[1]-k2c),2)}  k3 d={mp.nstr(abs(kk[2]-k3c),2)}")
    print(f"    k1={mp.nstr(k1c,14)} k2={mp.nstr(k2c,14)} k3={mp.nstr(k3c,14)}")
