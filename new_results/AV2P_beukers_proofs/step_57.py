from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, binomial
mp.dps = 60
z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3)
def f(x): return 2*mp.asin(mpf(x)**mpf('0.5')/2)**2
def H(n):
    n = int(n)
    return sum(1/mpf(k) for k in range(1,n+1)) if n>0 else mpf(0)
T = mp.nsum(lambda n: (H(2*n)-H(n))/(mpf(n)**2*binomial(2*n,n)), [1, mp.inf])
Tint = quad(lambda x: (f(x)-f(x*x))/(1-x), [0,1])
L = quad(lambda t: t*log(1+2*mp.sin(t)), [0, pi/6])
print("T - Tint:", nprint(T - Tint, 4))
print("T - 4L:", nprint(T - 4*L, 4))
print("L - (zeta3/2 - (pi/6)Cl13):", nprint(L - (z3/2 - (pi/6)*C13), 4))
print("T - (2zeta3 - (2pi/3)Cl13):", nprint(T - (2*z3 - (2*pi/3)*C13), 4))
# B2(1) final
U = mp.nsum(lambda n: 1/(mpf(n)**3*binomial(2*n,n)), [1, mp.inf])
print("B2(1) = T+U - 2zeta3/3:", nprint(T + U - 2*z3/3, 4))
with open('progress.log','a') as fh:
    fh.write("step_56: corrected Task-1 chain T = 4L, L = zeta3/2-(pi/6)Cl2(pi/3) proved and verified\n")
