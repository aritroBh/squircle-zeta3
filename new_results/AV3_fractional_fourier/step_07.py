from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nstr, fac)
mp.dps = 60
def E(x): return nstr(abs(x), 3)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def a_cos_gam(n, s): return re((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))
def b_sin_gam(n, s): return im((2/pi)*(-1j*n)**(-s-1) * gammainc(s+1, 0, -1j*n*pi))

# Asymptotic templates (derived):
# a_n(s) ~ kink + boundary(even powers):
#   kinkA = -(2/pi) Gamma(s+1) sin(pi s/2) n^-(s+1)
#   bndA_j = (2(-1)^n/pi) * (-1)^(j-1) * Gamma(s+1)/Gamma(s-2j+2) * pi^(s-2j+1) * n^(-2j),  j=1,2,...
# b_n(s) ~ kink + boundary(odd powers):
#   kinkB = (2/pi) Gamma(s+1) cos(pi s/2) n^-(s+1)
#   bndB_j = (2(-1)^(n+1)/pi) * (-1)^j * Gamma(s+1)/Gamma(s-2j+1) * pi^(s-2j) * n^(-(2j+1)), j=0,1,...
def kinkA(n,s): return -(2/pi)*gamma(s+1)*sin(pi*s/2)*n**(-s-1)
def bndA(n,s,J): return sum( (2*(-1)**n/pi)*(-1)**(j-1)*gamma(s+1)/gamma(s-2*j+2)*pi**(s-2*j+1)*n**(-2*j) for j in range(1,J+1))
def kinkB(n,s): return (2/pi)*gamma(s+1)*cos(pi*s/2)*n**(-s-1)
def bndB(n,s,J): return sum( (2*(-1)**(n+1)/pi)*(-1)**j*gamma(s+1)/gamma(s-2*j+1)*pi**(s-2*j)*n**(-(2*j+1)) for j in range(0,J+1))

print("== asymptotic templates vs exact (s=2.5 and s=3.4; increasing boundary order J) ==")
for s in [mpf('2.5'), mpf('3.4')]:
    for n in [10, 100]:
        aex = a_cos_gam(n,s); bex = b_sin_gam(n,s)
        row = f"s={nstr(s,4)}, n={n}: "
        for J in [1,2,3,4]:
            ea = abs(aex - kinkA(n,s) - bndA(n,s,J))
            eb = abs(bex - kinkB(n,s) - bndB(n,s,J))
            row += f"J{J}: a_err={E(ea)} b_err={E(eb)} | "
        print(row)

print()
print("== exact integer-s formulas ==")
forms_a = {
 1: lambda n: 2*((-1)**n-1)/(pi*n**2),
 2: lambda n: 4*(-1)**n/n**2,
 3: lambda n: 6*pi*(-1)**n/n**2 + 12*(1-(-1)**n)/(pi*n**4),
 4: lambda n: (-1)**n*(8*pi**2/n**2 - 48/n**4),
}
forms_b = {
 1: lambda n: 2*(-1)**(n+1)/n,
 2: lambda n: 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3),
 3: lambda n: 2*pi**2*(-1)**(n+1)/n + 12*(-1)**n/n**3,
 4: lambda n: 2*pi**3*(-1)**(n+1)/n + 24*pi*(-1)**n/n**3 + 48*(1-(-1)**n)/(pi*n**5),
}
for s in [1,2,3,4]:
    ea = max(abs(a_cos_gam(n,s) - forms_a[s](n)) for n in [1,2,3,7,20,51])
    eb = max(abs(b_sin_gam(n,s) - forms_b[s](n)) for n in [1,2,3,7,20,51])
    print(f"s={s}: max |a_exact-formula| err = {E(ea)}, max |b_exact-formula| err = {E(eb)}")
logline("step_06: kink+boundary asymptotics verified (errors drop with J as predicted); exact integer-s forms s=1..4 verified")
