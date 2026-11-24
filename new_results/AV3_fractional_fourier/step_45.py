from mpmath import (mp, mpf, pi, gamma, sin, cos, nstr, gammainc, re, im)
mp.dps = 80
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
def gsafe(z):
    try: return gamma(z)
    except ValueError: return mpf('inf')
def kinkA(n,s): return -(2/pi)*gamma(s+1)*sin(pi*s/2)*n**(-s-1)
def kinkB(n,s): return (2/pi)*gamma(s+1)*cos(pi*s/2)*n**(-s-1)
def bndA(n,s,J):
    tot = 0
    for j in range(1,J+1):
        g = gsafe(s-2*j+2)
        if g == mpf('inf'): continue
        tot += (2*(-1)**n/pi)*(-1)**(j-1)*gamma(s+1)/g*pi**(s-2*j+1)*n**(-2*j)
    return tot
def bndB(n,s,J):
    tot = 0
    for j in range(0,J+1):
        g = gsafe(s-2*j+1)
        if g == mpf('inf'): continue
        tot += (2*(-1)**(n+1)/pi)*(-1)**j*gamma(s+1)/g*pi**(s-2*j)*n**(-(2*j+1))
    return tot
def a_gam(n, s): return re((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))
def b_gam(n, s): return im((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))
nn = 300
print("s : a-kink pred vs n^{s+1}(a_n - bnd) ; b-kink pred vs got")
for s in [mpf(2), mpf('2.25'), mpf('2.5'), mpf(3), mpf('3.5'), mpf(4)]:
    ra = (a_gam(nn,s) - bndA(nn,s,2))*nn**(s+1)
    rb = (b_gam(nn,s) - bndB(nn,s,2))*nn**(s+1)
    print(f"s={nstr(s,5)}: a {nstr(kinkA(nn,s)*nn**(s+1),10)} | {nstr(ra,10)} ; b {nstr(kinkB(nn,s)*nn**(s+1),10)} | {nstr(rb,10)}")
logline("step_40: kink census with boundary subtraction, integer+half-integer s")
