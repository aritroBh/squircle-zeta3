import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, cos, sin, exp, matrix, lu_solve
mp.dps = 75

# --- cross-check P_squircle via two independent methods ---
def P_sq_support(tau):
    p = 1/(1-tau); q = 1/tau
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/4])
def P_sq_arclength(tau):
    p = 1/(1-tau)
    a = 2/p
    f = lambda t: sqrt( (a*cos(t)**(a-1)*sin(t))**2 + (a*sin(t)**(a-1)*cos(t))**2 )
    return 4*mpmath.quad(f, [0, pi/2])
t = mpf('0.02')
P1 = P_sq_support(t); P2 = P_sq_arclength(t)
print("support:", mp.nstr(P1,60)); print("arclen :", mp.nstr(P2,60))
print("diff:", mp.nstr(abs(P1-P2),5))

# --- robust coefficient extraction: interpolation on 6 well-separated points ---
taus6 = [mpf('0.032')/2**j for j in range(6)]   # 0.032 ... 0.001
Ps6 = [P_sq_support(t) for t in taus6]
P0 = 4*sqrt(2)
# divide out tau^2 first: Q(tau) = (P-P0)/tau^2 = c2 + c3 tau + ... + c7 tau^5
Q = [(P-P0)/t**2 for P,t in zip(Ps6,taus6)]
n=6
A = matrix(n,n)
for i,t in enumerate(taus6):
    for j in range(n):
        A[i,j] = t**j
c = lu_solve(A, matrix(n,1,lambda i,j: Q[i]))
names = ['c2','c3','c4','c5','c6','c7']
for nm,ci in zip(names,c):
    print(nm, "=", mp.nstr(ci,45))
print("pred c2 = sqrt2*zeta2 =", mp.nstr(sqrt(2)*zeta(2),45))
print("pred c3 = sqrt2*zeta3 =", mp.nstr(sqrt(2)*zeta(3),45))
